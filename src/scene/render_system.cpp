#include "scene/render_system.h"
#include "core/assert.h"
#include "core/bit.h"
#include "core/octahedral.h"
#include "core/thread_context.h"
#include "core/time.h"
#include "device/device.h"
#include "render/virtual_texture_manager.h"
#include "scene/clusterizer.h"
#include "scene/material.h"
#include "scene/scene.h"
#include "scene/shader_compilation.h"
#include <cstring>
#include <memory>
#include <tbb/parallel_for.h>
#include <unordered_set>

#define YBI_USE_CLUSTERS

namespace ybi
{

MaterialManager::~MaterialManager() = default;
RenderSystem::~RenderSystem() = default;

static bool MeshClusterDataEqual(const Mesh &a, const Mesh &b)
{
    if (a.GetPositionsHandle().index != b.GetPositionsHandle().index ||
        a.GetNormalsHandle().index != b.GetNormalsHandle().index ||
        a.GetUVsHandle().index != b.GetUVsHandle().index ||
        a.GetUVIndicesHandle().index != b.GetUVIndicesHandle().index ||
        a.GetIndicesHandle().index != b.GetIndicesHandle().index ||
        a.GetFaceMaterialIndicesHandle().index != b.GetFaceMaterialIndicesHandle().index ||
        a.GetMaterialIndex() != b.GetMaterialIndex() || a.GetUVType() != b.GetUVType() ||
        a.GetNormalType() != b.GetNormalType())
    {
        return false;
    }

    const std::vector<Primvar> &aPrimvars = a.GetPrimvars();
    const std::vector<Primvar> &bPrimvars = b.GetPrimvars();
    if (aPrimvars.size() != bPrimvars.size())
    {
        return false;
    }

    for (size_t i = 0; i < aPrimvars.size(); ++i)
    {
        const Primvar &aPrimvar = aPrimvars[i];
        const Primvar &bPrimvar = bPrimvars[i];
        if (aPrimvar.id != bPrimvar.id || aPrimvar.type != bPrimvar.type ||
            aPrimvar.valueType != bPrimvar.valueType || aPrimvar.index != bPrimvar.index ||
            aPrimvar.indicesIndex != bPrimvar.indicesIndex)
        {
            return false;
        }
    }
    return true;
}

static uint64_t HashMeshClusterData(const Mesh &mesh)
{
    uint64_t hash = 0;
    hash = HashCombine(hash, mesh.GetPositionsHandle().index);
    hash = HashCombine(hash, mesh.GetNormalsHandle().index);
    hash = HashCombine(hash, mesh.GetUVsHandle().index);
    hash = HashCombine(hash, mesh.GetUVIndicesHandle().index);
    hash = HashCombine(hash, mesh.GetIndicesHandle().index);
    hash = HashCombine(hash, mesh.GetFaceMaterialIndicesHandle().index);
    hash = HashCombine(hash, mesh.GetMaterialIndex());
    hash = HashCombine(hash, static_cast<uint64_t>(mesh.GetUVType()));
    hash = HashCombine(hash, static_cast<uint64_t>(mesh.GetNormalType()));
    for (const Primvar &primvar : mesh.GetPrimvars())
    {
        hash = HashCombine(hash, primvar.id);
        hash = HashCombine(hash, static_cast<uint64_t>(primvar.type));
        hash = HashCombine(hash, static_cast<uint64_t>(primvar.valueType));
        hash = HashCombine(hash, primvar.index);
        hash = HashCombine(hash, primvar.indicesIndex);
    }
    return hash;
}

RenderSystem::RenderSystem(Device *device,
                           DisplayDevice *displayDevice,
                           DataContainer *container,
                           Scene *rootScene,
                           VirtualTextureManager *virtualTextureManager)
    : device(device), displayDevice(displayDevice), dataContainer(container), rootScene(rootScene)
{
    deviceArena = device->CreateArena();
    materialManager.virtualTextureManager = virtualTextureManager;
}

void RenderSystem::Update(SceneImportContext *context)
{
    YBI_ASSERT(dataContainer);
    YBI_ASSERT(rootScene);

    materialManager.Update(device, displayDevice, dataContainer, &deviceScene, context);
    geometryManager.Update(device, deviceArena.get(), dataContainer, rootScene, &deviceScene);
}

constexpr size_t INVALID_INDEX = ~size_t(0);
constexpr uint32_t INVALID_OFFSET = ~0u;

template <typename T>
static size_t CalculateBufferOffsets(const std::vector<Array<T>> &buffers,
                                     std::vector<uint32_t> &offsets)
{
    offsets.resize(buffers.size(), INVALID_OFFSET);

    size_t totalSize = 0;
    for (size_t index = 0; index < buffers.size(); index++)
    {
        const Array<T> &buffer = buffers[index];
        if (buffer.size())
        {
            offsets[index] = static_cast<uint32_t>(totalSize);
            totalSize += buffer.size();
        }
    }

    return totalSize;
}

static size_t CalculateCompressedClusterByteOffsets(const std::vector<Array<uint8_t>> &buffers,
                                                    std::vector<uint32_t> &offsets)
{
    offsets.resize(buffers.size(), INVALID_OFFSET);

    size_t totalSize = 0;
    for (size_t index = 0; index < buffers.size(); index++)
    {
        const Array<uint8_t> &buffer = buffers[index];
        if (buffer.size())
        {
            totalSize = AlignUp(totalSize, sizeof(uint32_t));
            offsets[index] = static_cast<uint32_t>(totalSize);
            totalSize += buffer.size() + sizeof(uint32_t);
        }
    }

    return totalSize;
}

template <typename T>
static void UploadDeviceBuffers(Device *device,
                                T *&devicePtr,
                                const std::vector<Array<T>> &buffers,
                                const std::vector<uint32_t> &offsets,
                                size_t totalSize)
{
    devicePtr = nullptr;
    if (totalSize == 0)
    {
        return;
    }

    YBI_ASSERT(offsets.size() == buffers.size());

    Array<T> flattenedBuffers(totalSize);
    memset(flattenedBuffers.data(), 0, sizeof(T) * flattenedBuffers.size());
    tbb::parallel_for(size_t(0), buffers.size(), [&](size_t index) {
        const Array<T> &buffer = buffers[index];
        if (!buffer.size())
        {
            return;
        }

        YBI_ASSERT(offsets[index] != INVALID_OFFSET);
        YBI_ASSERT(static_cast<size_t>(offsets[index]) + buffer.size() <= totalSize);
        memcpy(flattenedBuffers.data() + offsets[index], buffer.data(), sizeof(T) * buffer.size());
    });

    devicePtr = reinterpret_cast<T *>(device->Alloc(sizeof(T) * totalSize).Data());
    device->CopyToDevice(devicePtr, flattenedBuffers.data(), sizeof(T) * flattenedBuffers.size());
}

static void UploadEncodedNormalBuffers(Device *device,
                                       uint32_t *&devicePtr,
                                       const std::vector<Array<Vec3f>> &buffers,
                                       const std::vector<uint32_t> &offsets,
                                       size_t totalSize)
{
    devicePtr = nullptr;
    if (totalSize == 0)
    {
        return;
    }

    Array<uint32_t> encodedNormals(totalSize);
    tbb::parallel_for(size_t(0), buffers.size(), [&](size_t index) {
        const Array<Vec3f> &buffer = buffers[index];
        if (!buffer.size())
        {
            return;
        }

        uint32_t *dst = encodedNormals.data() + offsets[index];
        for (size_t i = 0; i < buffer.size(); i++)
        {
            Vec3f v = buffer[i];
            // TODO: properly handle invalid normals
            if (LengthSquared(buffer[i]) < 1e-20f)
            {
                v = Vec3f(0.f, 0.f, 1.f);
            }
            dst[i] = EncodeOctahedral(Normalize(v));
        }
    });

    devicePtr = reinterpret_cast<uint32_t *>(device->Alloc(sizeof(uint32_t) * totalSize).Data());
    device->CopyToDevice(devicePtr, encodedNormals.data(), sizeof(uint32_t) * totalSize);
}

static void
BuildTLAS(Device *device, DeviceArena *deviceArena, DataContainer *container, Scene *rootScene)
{
    YBI_ASSERT(device);
    YBI_ASSERT(deviceArena);
    YBI_ASSERT(container);
    YBI_ASSERT(rootScene);

    // TODO: supports 2+ level scenes
    container->GetBVHs(GeometryType::Scene)
        .EmplaceBack(device->BuildInstancesBVH(deviceArena, rootScene->instances, container));
}

// GeometryManager
void GeometryManager::Update(Device *device,
                             DeviceArena *deviceArena,
                             DataContainer *container,
                             Scene *rootScene,
                             DeviceScene *deviceScene)
{
    // TODO: this function is a bit overloaded maybe
    ScratchArena scratch;

    // Set up device scene
    View<Object> objects = scratch.arena->PushArray<Object>(container->meshes.size());
    View<DeviceMaterial> hostMaterials;
    View<MaterialHandle> hostMaterialHandles;
    if (container->materials.size())
    {
        hostMaterials = scratch.arena->PushArray<DeviceMaterial>(container->materials.size());
    }
    if (container->materialHandles.size())
    {
        hostMaterialHandles =
            scratch.arena->PushArray<MaterialHandle>(container->materialHandles.size());
        memcpy(hostMaterialHandles.data(),
               container->materialHandles.data(),
               sizeof(container->materialHandles[0]) * container->materialHandles.size());
    }

    // Prepare device materials
    for (size_t materialIndex = 0; materialIndex < container->materials.size(); materialIndex++)
    {
        const Material &material = container->materials[materialIndex];
        DeviceMaterial &deviceMaterial = hostMaterials[materialIndex];
        deviceMaterial.udimBase = material.udimBase;
        deviceMaterial.udimTableOffset = material.udimTableOffset;
        deviceMaterial.valid = material.valid;
    }

    std::vector<uint32_t> positionOffsets;
    std::vector<uint32_t> normalOffsets;
    std::vector<uint32_t> uvOffsets;
    std::vector<uint32_t> uvIndexOffsets;
    std::vector<uint32_t> triangleIndexOffsets;
    std::vector<uint32_t> triangleMaterialIndexOffsets;
    std::vector<uint32_t> meshPrimvarOffsets;

    size_t totalNumVertices = CalculateBufferOffsets(container->meshPositions, positionOffsets);
    size_t totalNumNormals = CalculateBufferOffsets(container->meshNormals, normalOffsets);
    size_t totalNumUVs = CalculateBufferOffsets(container->meshUVs, uvOffsets);
    size_t totalNumUVIndices = CalculateBufferOffsets(container->meshUVIndices, uvIndexOffsets);
    size_t totalNumTriangleIndices =
        CalculateBufferOffsets(container->meshTriangleIndices, triangleIndexOffsets);
    size_t totalNumTriangleMaterialIndices =
        CalculateBufferOffsets(container->meshFaceMaterialIndices, triangleMaterialIndexOffsets);
    meshPrimvarOffsets.resize(container->meshes.size());

    // Calculate primvar offsets
    size_t primvarOffset = 0;
    for (size_t meshIndex = 0; meshIndex < container->meshes.size(); meshIndex++)
    {
        const Mesh &mesh = container->meshes[meshIndex];
        size_t numPrimvars = mesh.GetPrimvars().size();

        meshPrimvarOffsets[meshIndex] = primvarOffset;
        primvarOffset += numPrimvars;
    }

    std::vector<size_t> bvhSourceIndices;
    bvhSourceIndices.reserve(container->meshes.size());
    std::unordered_map<uint64_t, size_t> bvhPairToMeshIndex;

    // Unlike a BVH, a compressed cluster embeds shading data. Reuse it only when all
    // standard mesh buffers, material assignments, and primvars match.
    std::vector<size_t> clusterSourceIndices;
    clusterSourceIndices.reserve(container->meshes.size());
    std::unordered_multimap<uint64_t, size_t> clusterHashToMeshIndex;

    for (size_t meshIndex = 0; meshIndex < container->meshes.size(); meshIndex++)
    {
        Mesh &mesh = container->meshes[meshIndex];

        Object &object = objects[meshIndex];
        object.positionOffset = INVALID_OFFSET;
        object.normalOffset = INVALID_OFFSET;
        object.uvOffset = INVALID_OFFSET;
        object.uvIndexOffset = INVALID_OFFSET;
        object.triangleIndexOffset = INVALID_OFFSET;
        object.triangleMaterialOffset = INVALID_OFFSET;
        object.clusterOffset = INVALID_OFFSET;
        object.uvType = mesh.GetUVType();
        object.normalType = mesh.GetNormalType();
        object.materialIndex = -1;
        object.primvarOffset = meshPrimvarOffsets[meshIndex];
        object.primvarCount = mesh.GetPrimvars().size();

        BufferHandle positionHandle = mesh.GetPositionsHandle();
        BufferHandle normalHandle = mesh.GetNormalsHandle();
        BufferHandle uvHandle = mesh.GetUVsHandle();
        BufferHandle uvIndexHandle = mesh.GetUVIndicesHandle();
        BufferHandle indexHandle = mesh.GetIndicesHandle();
        BufferHandle triangleMaterialHandle = mesh.GetFaceMaterialIndicesHandle();

        if (positionHandle.IsValid())
        {
            object.positionOffset = positionOffsets[positionHandle.index];
        }

        if (normalHandle.IsValid())
        {
            object.normalOffset = normalOffsets[normalHandle.index];
        }

        if (indexHandle.IsValid())
        {
            object.triangleIndexOffset = triangleIndexOffsets[indexHandle.index];
        }

        if (uvHandle.IsValid())
        {
            object.uvOffset = uvOffsets[uvHandle.index];
        }

        if (uvIndexHandle.IsValid())
        {
            object.uvIndexOffset = uvIndexOffsets[uvIndexHandle.index];
        }

        if (triangleMaterialHandle.IsValid())
        {
            object.triangleMaterialOffset =
                triangleMaterialIndexOffsets[triangleMaterialHandle.index];
        }

        uint32_t materialIndex = mesh.GetMaterialIndex();
        if (materialIndex < container->materialHandles.size())
        {
            object.materialIndex = static_cast<int>(materialIndex);
        }

        if (positionHandle.IsValid() && indexHandle.IsValid())
        {
            uint64_t key = (uint64_t(positionHandle.index) << 32) | uint64_t(indexHandle.index);
            auto it = bvhPairToMeshIndex.find(key);
            if (it == bvhPairToMeshIndex.end())
            {
                bvhPairToMeshIndex.insert({key, meshIndex});
                bvhSourceIndices.push_back(meshIndex);
            }
            else
            {
                bvhSourceIndices.push_back(it->second);
            }
        }
        else
        {
            bvhSourceIndices.push_back(meshIndex);
        }

        const uint64_t clusterHash = HashMeshClusterData(mesh);
        const auto clusterCandidates = clusterHashToMeshIndex.equal_range(clusterHash);
        size_t clusterSourceIndex = meshIndex;
        for (auto it = clusterCandidates.first; it != clusterCandidates.second; ++it)
        {
            if (MeshClusterDataEqual(mesh, container->meshes[it->second]))
            {
                clusterSourceIndex = it->second;
                break;
            }
        }
        if (clusterSourceIndex == meshIndex)
        {
            clusterHashToMeshIndex.emplace(clusterHash, meshIndex);
        }
        clusterSourceIndices.push_back(clusterSourceIndex);
    }

    // Create BVHs
    container->GetBVHs(GeometryType::Mesh).Reserve(container->meshes.size());
    printf("num meshes: %zu\n", container->meshes.size());

#ifdef YBI_USE_CLUSTERS
    if (device->SupportsClusterBVHs())
    {
        std::vector<MeshletClustersResult> clusterResults(container->meshes.size());
        std::vector<Array<DeviceCluster>> compressedMeshClusters(container->meshes.size());

        tbb::parallel_for(size_t(0), container->meshes.size(), [&](size_t meshIndex) {
            Mesh &mesh = container->meshes[meshIndex];
            const size_t sourceIndex = clusterSourceIndices[meshIndex];

            if (sourceIndex != meshIndex)
            {
                return;
            }

            if (!mesh.GetPositionsHandle().IsValid() || !mesh.GetIndicesHandle().IsValid())
            {
                return;
            }

            const Array<int> &sourceIndices =
                container->meshTriangleIndices[mesh.GetIndicesHandle().index];
            YBI_ASSERT(sourceIndices.size() % 3 == 0);

            Array<int> filteredIndices;
            filteredIndices.Reserve(sourceIndices.size());
            for (size_t i = 0; i < sourceIndices.size(); i += 3)
            {
                const int index0 = sourceIndices[i + 0];
                const int index1 = sourceIndices[i + 1];
                const int index2 = sourceIndices[i + 2];
                if (index0 != index1 && index1 != index2 && index2 != index0)
                {
                    filteredIndices.EmplaceBack(index0);
                    filteredIndices.EmplaceBack(index1);
                    filteredIndices.EmplaceBack(index2);
                }
            }

            MeshletClustersResult &result = clusterResults[meshIndex];
            ClusterizeMesh(container, mesh, filteredIndices, result);
            compressedMeshClusters[meshIndex].Resize(result.meshlets.size());
        });

        std::vector<uint32_t> clusterOffsets;
        const size_t totalNumClusters =
            CalculateBufferOffsets(compressedMeshClusters, clusterOffsets);

        container->compressedClusterBytes.clear();
        container->compressedClusterBytes.resize(totalNumClusters);
        const Timer compressClustersTimer = StartTimer();
        tbb::parallel_for(size_t(0), container->meshes.size(), [&](size_t meshIndex) {
            if (clusterSourceIndices[meshIndex] != meshIndex)
            {
                return;
            }

            YBI_ASSERT(clusterOffsets[meshIndex] != INVALID_OFFSET);
            CompressClusters(clusterResults[meshIndex],
                             container->meshes[meshIndex],
                             container,
                             compressedMeshClusters[meshIndex],
                             clusterOffsets[meshIndex]);
        });
        const double compressClustersMilliseconds = StopTimer(compressClustersTimer);
        printf("CompressClusters parallel_for: %.3f ms\n", compressClustersMilliseconds);

        std::vector<uint32_t> compressedClusterByteOffsets;
        const size_t totalCompressedClusterBytes = CalculateCompressedClusterByteOffsets(
            container->compressedClusterBytes, compressedClusterByteOffsets);

        for (size_t meshIndex = 0; meshIndex < container->meshes.size(); meshIndex++)
        {
            if (clusterSourceIndices[meshIndex] != meshIndex)
            {
                continue;
            }

            YBI_ASSERT(clusterOffsets[meshIndex] != INVALID_OFFSET);
            Array<DeviceCluster> &clusters = compressedMeshClusters[meshIndex];
            for (size_t clusterIndex = 0; clusterIndex < clusters.size(); clusterIndex++)
            {
                const size_t globalClusterIndex = clusterOffsets[meshIndex] + clusterIndex;
                YBI_ASSERT(globalClusterIndex < compressedClusterByteOffsets.size());
                clusters[clusterIndex].byteOffset =
                    compressedClusterByteOffsets[globalClusterIndex];
            }
        }

        for (size_t meshIndex = 0; meshIndex < container->meshes.size(); meshIndex++)
        {
            const size_t sourceIndex = clusterSourceIndices[meshIndex];
            YBI_ASSERT(sourceIndex < clusterOffsets.size());
            YBI_ASSERT(clusterOffsets[sourceIndex] != INVALID_OFFSET);
            objects[meshIndex].triangleIndexOffset = INVALID_OFFSET;
            objects[meshIndex].clusterOffset = clusterOffsets[sourceIndex];
        }

        Array<BVH> &bvhs = container->GetBVHs(GeometryType::Mesh);
        const Timer buildClusterBVHsTimer = StartTimer();
        for (size_t meshIndex = 0; meshIndex < container->meshes.size(); meshIndex++)
        {
            Mesh &mesh = container->meshes[meshIndex];
            const size_t sourceIndex = bvhSourceIndices[meshIndex];

            if (sourceIndex != meshIndex)
            {
                bvhs.EmplaceBack(bvhs[sourceIndex]);
                continue;
            }

            const MeshletClustersResult &result = clusterResults[meshIndex];
            YBI_ASSERT(result.meshlets.size());
            bvhs.EmplaceBack(device->BuildClusterBVH(deviceArena, container, mesh, result));
        }
        const double buildClusterBVHsMilliseconds = StopTimer(buildClusterBVHsTimer);
        printf("BuildClusterBVHs: %.3f ms\n", buildClusterBVHsMilliseconds);

        BuildTLAS(device, deviceArena, container, rootScene);

        deviceScene->positions = nullptr;
        deviceScene->normals = nullptr;
        deviceScene->uvs = nullptr;
        deviceScene->uvIndices = nullptr;
        deviceScene->triangleIndices = nullptr;
        deviceScene->triangleMaterialIndices = nullptr;
        UploadDeviceBuffers(device,
                            deviceScene->compressedClusters,
                            compressedMeshClusters,
                            clusterOffsets,
                            totalNumClusters);
        UploadDeviceBuffers(device,
                            deviceScene->compressedClusterBytes,
                            container->compressedClusterBytes,
                            compressedClusterByteOffsets,
                            totalCompressedClusterBytes);
    }
    else
#endif
    {
        const Timer buildTriangleBVHsTimer = StartTimer();
        for (size_t meshIndex = 0; meshIndex < container->meshes.size(); meshIndex++)
        {
            Mesh &mesh = container->meshes[meshIndex];
            const size_t sourceIndex = bvhSourceIndices[meshIndex];

            if (sourceIndex != meshIndex)
            {
                Array<BVH> &bvhs = container->GetBVHs(GeometryType::Mesh);
                // TODO: this is brittle. also bvhs probably shouldn't be duplicated
                bvhs.EmplaceBack(bvhs[sourceIndex]);
            }
            else
            {
                container->GetBVHs(GeometryType::Mesh)
                    .EmplaceBack(device->BuildTriangleBVH(deviceArena, container, mesh));
            }
        }
        const double buildTriangleBVHsMilliseconds = StopTimer(buildTriangleBVHsTimer);
        printf("BuildTriangleBVHs: %.3f ms\n", buildTriangleBVHsMilliseconds);

        BuildTLAS(device, deviceArena, container, rootScene);

        UploadDeviceBuffers(device,
                            deviceScene->triangleIndices,
                            container->meshTriangleIndices,
                            triangleIndexOffsets,
                            totalNumTriangleIndices);
        UploadDeviceBuffers(device,
                            deviceScene->positions,
                            container->meshPositions,
                            positionOffsets,
                            totalNumVertices);
        UploadEncodedNormalBuffers(
            device, deviceScene->normals, container->meshNormals, normalOffsets, totalNumNormals);
        UploadDeviceBuffers(device, deviceScene->uvs, container->meshUVs, uvOffsets, totalNumUVs);
        UploadDeviceBuffers(device,
                            deviceScene->uvIndices,
                            container->meshUVIndices,
                            uvIndexOffsets,
                            totalNumUVIndices);
        UploadDeviceBuffers(device,
                            deviceScene->triangleMaterialIndices,
                            container->meshFaceMaterialIndices,
                            triangleMaterialIndexOffsets,
                            totalNumTriangleMaterialIndices);
        deviceScene->compressedClusters = nullptr;
        deviceScene->compressedClusterBytes = nullptr;
    }

    // Upload primvars
    if (primvarOffset)
    {
        std::vector<DeviceAttribute> attributes;
        for (size_t i = 0; i < container->meshes.size(); i++)
        {
            const Mesh &mesh = container->meshes[i];
            for (const Primvar &primvar : mesh.GetPrimvars())
            {
                DeviceAttribute deviceAttribute;
                deviceAttribute.id = primvar.id;
                deviceAttribute.attributeType = primvar.type;

                // TODO: fill these out when clusterizer is not used
                // deviceAttribute.indices = ;
                // deviceAttribute.indicesIndex = ;
                attributes.push_back(deviceAttribute);
            }
        }

        deviceScene->deviceAttributes = reinterpret_cast<DeviceAttribute *>(
            device->Alloc(sizeof(DeviceAttribute) * primvarOffset).Data());
        device->CopyToDevice(deviceScene->deviceAttributes,
                             attributes.data(),
                             sizeof(DeviceAttribute) * primvarOffset);
    }

    Help();

    printf("total num vertices: %zu\n", totalNumVertices);
    printf("num instances: %zu\n", rootScene->instances.size());

    deviceScene->materials = nullptr;
    if (hostMaterials.size())
    {
        deviceScene->materials = reinterpret_cast<DeviceMaterial *>(
            device->Alloc(sizeof(DeviceMaterial) * hostMaterials.size()).Data());
        device->CopyToDevice(
            deviceScene->materials, hostMaterials.data(), hostMaterials.numBytes());
    }
    if (hostMaterialHandles.size())
    {
        deviceScene->materialHandles = reinterpret_cast<MaterialHandle *>(
            device->Alloc(sizeof(MaterialHandle) * hostMaterialHandles.size()).Data());
        device->CopyToDevice(deviceScene->materialHandles,
                             hostMaterialHandles.data(),
                             hostMaterialHandles.numBytes());
    }

    deviceScene->objects = reinterpret_cast<Object *>(
        device->Alloc(sizeof(Object) * container->meshes.size()).Data());
    device->CopyToDevice(deviceScene->objects, objects.data(), objects.numBytes());

    printf("bvh total allocated: %zu\n", device->GetBVHAllocationSize());
}

void MaterialManager::AddTextures(OpenPbrMaterial &material, MaterialNodeId member)
{
    const MaterialNode *node = &material.graph.nodes[member];

    // TODO: allow enable/disable vt system
    if (IsTextureOpcode(node->type))
    {
        bool found = false;
        for (const MaterialInput &input : node->inputs)
        {
            std::visit(
                [&](const auto &value) {
                    using T = std::decay_t<decltype(value)>;
                    if constexpr (std::is_same_v<T, UdimTextureVector>)
                    {
                        MaterialTextureHandle handle = virtualTextureManager->AddTexture(value);
                        found = true;
                        textureNodeToHandle[node] = handle;
                    }
                },
                input.value);
            if (found)
                break;
        }
        YBI_ASSERT(found);
    }
}

DeviceMaterialExpression MaterialManager::ProcessMaterialGraph(OpenPbrMaterial &material,
                                                               const char *memberName,
                                                               MaterialNodeId member,
                                                               std::vector<uint32_t> &nodes)
{
    (void)memberName;
    return CompileShader(
        material.graph, member, textureNodeToHandle, attributeManager.attributeIdMap, nodes);
}

AttributeId AttributeManager::GetOrAdd(const std::string &name)
{
    auto find = attributeIdMap.find(name);
    if (find != attributeIdMap.end())
    {
        return find->second;
    }
    else
    {
        AttributeId id = counter++;
        attributeIdMap[name] = id;
        return id;
    }
}

void AttributeManager::FindPrimvars(DataContainer *container, SceneImportContext *context)
{
    std::vector<std::unordered_set<std::string>> materialGeomProps;

    for (OpenPbrMaterial &material : container->openPbrMaterials)
    {
        materialGeomProps.emplace_back();
        std::unordered_set<std::string> &geomProps = materialGeomProps.back();

        for (MaterialNode &node : material.graph.nodes)
        {
            if (IsAttributeOpcode(node.type))
            {
                for (const MaterialInput &input : node.inputs)
                {
                    if (input.name == "geomprop")
                    {
                        const std::string *value = std::get_if<std::string>(&input.value);
                        YBI_ASSERT(value != nullptr);

                        geomProps.insert(*value);
                    }
                }
            }
        }
    }

    for (PendingMesh &pendingMesh : context->pendingMeshes)
    {
        Mesh &mesh = container->meshes[pendingMesh.meshIndex];
        YBI_ASSERT(mesh.HasMaterial());

        std::vector<int> materialIndices;

        if (mesh.HasSingleMaterial())
        {
            materialIndices.push_back(mesh.GetMaterialIndex());
        }
        else if (mesh.GetFaceMaterialIndicesHandle().IsValid())
        {
            std::unordered_set<int> uniqueMaterialIndexSet;
            const Array<int> &indices =
                container->meshFaceMaterialIndices[mesh.GetFaceMaterialIndicesHandle().index];
            for (int index : indices)
            {
                if (index >= 0)
                {
                    uniqueMaterialIndexSet.insert(index);
                }
            }
            materialIndices.insert(materialIndices.end(),
                                   uniqueMaterialIndexSet.begin(),
                                   uniqueMaterialIndexSet.end());
        }

        for (int materialIndex : materialIndices)
        {
            const MaterialHandle &handle = container->materialHandles[materialIndex];
            if (!handle.isGraph)
            {
                continue;
            }

            const std::unordered_set<std::string> &geomProps = materialGeomProps[handle.index];
            for (const std::string &name : geomProps)
            {
                const AttributeId id = GetOrAdd(name);

                // TODO: these asserts should probably be handled as edge cases instead of
                // crashing
                YBI_ASSERT(mesh.GetPrimvar(id) == nullptr);

                MeshAttribute attribute;
                bool result = pendingMesh.loadAttribute(name, attribute);

                YBI_ASSERT(result);

                Primvar primvar = container->attributes.AddAttribute(
                    id, attribute.type, std::move(attribute.value), std::move(attribute.indices));
                mesh.SetPrimvar(std::move(primvar));
            }
        }
    }

    context->pendingMeshes.clear();
}

void MaterialManager::Update(Device *device,
                             DisplayDevice *displayDevice,
                             DataContainer *container,
                             DeviceScene *deviceScene,
                             SceneImportContext *importContext)
{
    // what's next
    // 1. uploading primvars to gpu (accounting for clusterizer...) <-- fin
    // 2. uploading mesh primvar handles to gpu (so shader can lookup right id) <-- fin
    // 3. runtime access of primvars <-- fin
    // 4. runtime shading

    attributeManager.FindPrimvars(container, importContext);

    YBI_ASSERT(virtualTextureManager != nullptr);

    // Convert textures to handles
    for (OpenPbrMaterial &material : container->openPbrMaterials)
    {
        for (MaterialNodeId nodeId = 0;
             nodeId < static_cast<MaterialNodeId>(material.graph.nodes.size());
             nodeId++)
        {
            AddTextures(material, nodeId);
        }
    }

    // Compile shader nodes
    std::vector<uint32_t> shaderNodes;

    ScratchArena scratch;
    View<DeviceMaterial2> hostMaterials;
    if (container->openPbrMaterials.size())
    {
        hostMaterials =
            scratch.arena->PushArray<DeviceMaterial2>(container->openPbrMaterials.size());
    }

    for (size_t materialIndex = 0; materialIndex < container->openPbrMaterials.size();
         materialIndex++)
    {
        OpenPbrMaterial &material = container->openPbrMaterials[materialIndex];
        DeviceMaterial2 &deviceMaterial = hostMaterials[materialIndex];
        deviceMaterial.base_color =
            ProcessMaterialGraph(material, "base_color", material.base_color, shaderNodes);
        // #define VISIT_PARAMETER(name, type) \
        //     deviceMaterial.name = ProcessMaterialGraph(material, #name, material.name,
        //     shaderNodes);
        //
        //         OPENPBR_PARAMETERS(VISIT_PARAMETER)
        // #undef VISIT_PARAMETER
    }

    if (hostMaterials.size())
    {
        deviceScene->deviceMaterial2s = reinterpret_cast<DeviceMaterial2 *>(
            device->Alloc(sizeof(DeviceMaterial2) * hostMaterials.size()).Data());
        device->CopyToDevice(
            deviceScene->deviceMaterial2s, hostMaterials.data(), hostMaterials.numBytes());
    }

    deviceScene->shaderNodes = nullptr;
    if (!shaderNodes.empty())
    {
        deviceScene->shaderNodes = reinterpret_cast<uint32_t *>(
            device->Alloc(sizeof(uint32_t) * shaderNodes.size()).Data());
        device->CopyToDevice(
            deviceScene->shaderNodes, shaderNodes.data(), sizeof(uint32_t) * shaderNodes.size());
    }
}

} // namespace ybi
