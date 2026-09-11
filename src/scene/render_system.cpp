#include "scene/render_system.h"
#include "core/assert.h"
#include "core/bit.h"
#include "core/octahedral.h"
#include "core/thread_context.h"
#include "core/time.h"
#include "device/device.h"
#include "render/virtual_texture_manager.h"
#include "scene/clusterizer.h"
#include "scene/device_cluster.h"
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
    size_t aIndex = 0;
    size_t bIndex = 0;
    for (;;)
    {
        while (aIndex < aPrimvars.size() && !IsPrimvarArray(aPrimvars[aIndex].valueType))
        {
            ++aIndex;
        }
        while (bIndex < bPrimvars.size() && !IsPrimvarArray(bPrimvars[bIndex].valueType))
        {
            ++bIndex;
        }

        if (aIndex == aPrimvars.size() || bIndex == bPrimvars.size())
        {
            return aIndex == aPrimvars.size() && bIndex == bPrimvars.size();
        }

        const Primvar &aPrimvar = aPrimvars[aIndex++];
        const Primvar &bPrimvar = bPrimvars[bIndex++];
        if (aPrimvar.id != bPrimvar.id || aPrimvar.type != bPrimvar.type ||
            aPrimvar.valueType != bPrimvar.valueType || aPrimvar.index != bPrimvar.index ||
            aPrimvar.indicesIndex != bPrimvar.indicesIndex)
        {
            return false;
        }
    }
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
        if (!IsPrimvarArray(primvar.valueType))
        {
            continue;
        }

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

template <typename T>
static size_t CalculateBufferOffsetsForMeshes(const std::vector<Array<T>> &buffers,
                                              std::vector<uint32_t> &offsets,
                                              const std::vector<Mesh> &meshes,
                                              const Array<BufferHandle> &meshHandles,
                                              BufferHandle (Mesh::*getBufferHandle)() const)
{
    offsets.assign(buffers.size(), INVALID_OFFSET);

    size_t totalSize = 0;
    for (size_t i = 0; i < meshHandles.size(); i++)
    {
        const BufferHandle meshHandle = meshHandles[i];
        YBI_ASSERT(meshHandle.IsValid());
        YBI_ASSERT(meshHandle.index < meshes.size());

        const BufferHandle handle = (meshes[meshHandle.index].*getBufferHandle)();
        if (!handle.IsValid())
        {
            continue;
        }

        YBI_ASSERT(handle.index < buffers.size());

        if (offsets[handle.index] != INVALID_OFFSET)
        {
            continue;
        }

        const Array<T> &buffer = buffers[handle.index];
        if (buffer.size())
        {
            offsets[handle.index] = static_cast<uint32_t>(totalSize);
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
        if (!buffer.size() || offsets[index] == INVALID_OFFSET)
        {
            return;
        }

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

    YBI_ASSERT(offsets.size() == buffers.size());

    Array<uint32_t> encodedNormals(totalSize);
    tbb::parallel_for(size_t(0), buffers.size(), [&](size_t index) {
        const Array<Vec3f> &buffer = buffers[index];
        if (!buffer.size() || offsets[index] == INVALID_OFFSET)
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

static void BuildTLAS(Device *device,
                      DeviceArena *deviceArena,
                      DataContainer *container,
                      const Scene *rootScene,
                      const Array<int> &objectToBVHIndex)
{
    YBI_ASSERT(device);
    YBI_ASSERT(deviceArena);
    YBI_ASSERT(container);
    YBI_ASSERT(rootScene);

    // TODO: supports 2+ level scenes
    container->GetBVHs(GeometryType::Scene)
        .EmplaceBack(device->BuildInstancesBVH(
            deviceArena, rootScene->instances, container, objectToBVHIndex));
}

void GeometryManager::BuildTriangleMeshBVHsAndUpload(Device *device,
                                                     DeviceArena *deviceArena,
                                                     DataContainer *container,
                                                     DeviceScene *deviceScene,
                                                     const Array<BufferHandle> &meshHandles,
                                                     Array<Object> &objects,
                                                     Array<int> &objectBVHIndices)
{
    YBI_ASSERT(objects.size() == objectBVHIndices.size());
    std::vector<uint32_t> positionOffsets;
    std::vector<uint32_t> normalOffsets;
    std::vector<uint32_t> uvOffsets;
    std::vector<uint32_t> uvIndexOffsets;
    std::vector<uint32_t> triangleIndexOffsets;
    std::vector<uint32_t> triangleMaterialIndexOffsets;

    size_t totalNumVertices = CalculateBufferOffsetsForMeshes(container->meshPositions,
                                                              positionOffsets,
                                                              container->meshes,
                                                              meshHandles,
                                                              &Mesh::GetPositionsHandle);

    size_t totalNumNormals = CalculateBufferOffsetsForMeshes(container->meshNormals,
                                                             normalOffsets,
                                                             container->meshes,
                                                             meshHandles,
                                                             &Mesh::GetNormalsHandle);

    size_t totalNumUVs = CalculateBufferOffsetsForMeshes(
        container->meshUVs, uvOffsets, container->meshes, meshHandles, &Mesh::GetUVsHandle);

    size_t totalNumUVIndices = CalculateBufferOffsetsForMeshes(container->meshUVIndices,
                                                               uvIndexOffsets,
                                                               container->meshes,
                                                               meshHandles,
                                                               &Mesh::GetUVIndicesHandle);

    size_t totalNumTriangleIndices =
        CalculateBufferOffsetsForMeshes(container->meshTriangleIndices,
                                        triangleIndexOffsets,
                                        container->meshes,
                                        meshHandles,
                                        &Mesh::GetIndicesHandle);

    size_t totalNumTriangleMaterialIndices =
        CalculateBufferOffsetsForMeshes(container->meshFaceMaterialIndices,
                                        triangleMaterialIndexOffsets,
                                        container->meshes,
                                        meshHandles,
                                        &Mesh::GetFaceMaterialIndicesHandle);

    std::vector<size_t> bvhSourceIndices;
    bvhSourceIndices.reserve(meshHandles.size());
    std::unordered_map<uint64_t, size_t> bvhPairToHandleIndex;

    for (size_t handleIndex = 0; handleIndex < meshHandles.size(); handleIndex++)
    {
        const size_t meshIndex = meshHandles[handleIndex].index;
        Mesh &mesh = container->meshes[meshIndex];
        Object &object = objects[meshIndex];
        object.positionOffset = INVALID_OFFSET;
        object.normalOffset = INVALID_OFFSET;
        object.uvOffset = INVALID_OFFSET;
        object.uvIndexOffset = INVALID_OFFSET;
        object.triangleIndexOffset = INVALID_OFFSET;
        object.triangleMaterialOffset = INVALID_OFFSET;
        object.clusterOffset = INVALID_OFFSET;

        // TODO: support primvars in the non clustered path
        YBI_ASSERT(mesh.GetPrimvars().size() == 0);

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

        if (positionHandle.IsValid() && indexHandle.IsValid())
        {
            uint64_t key = (uint64_t(positionHandle.index) << 32) | uint64_t(indexHandle.index);
            auto it = bvhPairToHandleIndex.find(key);
            if (it == bvhPairToHandleIndex.end())
            {
                bvhPairToHandleIndex.insert({key, handleIndex});
                bvhSourceIndices.push_back(handleIndex);
            }
            else
            {
                bvhSourceIndices.push_back(it->second);
            }
        }
        else
        {
            bvhSourceIndices.push_back(handleIndex);
        }
    }

    Array<BVH> &bvhs = container->GetBVHs(GeometryType::Mesh);
    const Timer buildTriangleBVHsTimer = StartTimer();
    for (size_t handleIndex = 0; handleIndex < meshHandles.size(); handleIndex++)
    {
        const size_t meshIndex = meshHandles[handleIndex].index;
        Mesh &mesh = container->meshes[meshIndex];
        const size_t sourceIndex = bvhSourceIndices[handleIndex];

        if (sourceIndex != handleIndex)
        {
            const size_t sourceMeshIndex = meshHandles[sourceIndex].index;
            YBI_ASSERT(objectBVHIndices[sourceMeshIndex] >= 0);
            objectBVHIndices[meshIndex] = objectBVHIndices[sourceMeshIndex];
        }
        else
        {
            bvhs.EmplaceBack(device->BuildTriangleBVH(deviceArena, container, mesh));
            objectBVHIndices[meshIndex] = static_cast<int>(bvhs.size() - 1);
        }
    }

    YBI_ASSERT(objects.size() == objectBVHIndices.size());
    const double buildTriangleBVHsMilliseconds = StopTimer(buildTriangleBVHsTimer);
    printf("BuildTriangleBVHs: %.3f ms\n", buildTriangleBVHsMilliseconds);

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
}

void GeometryManager::BuildClusterMeshBVHsAndUpload(Device *device,
                                                    DeviceArena *deviceArena,
                                                    DataContainer *container,
                                                    DeviceScene *deviceScene,
                                                    const Array<BufferHandle> &meshHandles,
                                                    Array<Object> &objects,
                                                    Array<int> &objectBVHIndices)
{
    YBI_ASSERT(objects.size() == objectBVHIndices.size());
    // Contains index into meshHandles
    std::vector<size_t> clusterSourceIndices;
    std::unordered_multimap<uint64_t, size_t> clusterHashToHandleIndex;
    clusterSourceIndices.reserve(meshHandles.size());

    std::vector<size_t> bvhSourceIndices;
    bvhSourceIndices.reserve(meshHandles.size());
    std::unordered_map<uint64_t, size_t> bvhPairToHandleIndex;

    // If meshes have same position/indices
    // If meshes have exact same data, then reuse their cluster data
    for (size_t handleIndex = 0; handleIndex < meshHandles.size(); handleIndex++)
    {
        BufferHandle handle = meshHandles[handleIndex];
        const Mesh &mesh = container->meshes[handle.index];

        BufferHandle positionHandle = mesh.GetPositionsHandle();
        BufferHandle indexHandle = mesh.GetIndicesHandle();
        if (positionHandle.IsValid() && indexHandle.IsValid())
        {
            uint64_t key = (uint64_t(positionHandle.index) << 32) | uint64_t(indexHandle.index);
            auto it = bvhPairToHandleIndex.find(key);
            if (it == bvhPairToHandleIndex.end())
            {
                bvhPairToHandleIndex.insert({key, handleIndex});
                bvhSourceIndices.push_back(handleIndex);
            }
            else
            {
                bvhSourceIndices.push_back(it->second);
            }
        }
        else
        {
            bvhSourceIndices.push_back(handleIndex);
        }

        const uint64_t clusterHash = HashMeshClusterData(mesh);
        const auto clusterCandidates = clusterHashToHandleIndex.equal_range(clusterHash);
        size_t clusterSourceIndex = handleIndex;
        for (auto it = clusterCandidates.first; it != clusterCandidates.second; ++it)
        {
            const size_t candidateHandleIndex = it->second;
            const size_t candidateMeshIndex = meshHandles[candidateHandleIndex].index;

            if (MeshClusterDataEqual(mesh, container->meshes[candidateMeshIndex]))
            {
                clusterSourceIndex = candidateHandleIndex;
                break;
            }
        }
        if (clusterSourceIndex == handleIndex)
        {
            clusterHashToHandleIndex.emplace(clusterHash, handleIndex);
        }
        clusterSourceIndices.push_back(clusterSourceIndex);
    }

    std::vector<MeshletClustersResult> clusterResults(meshHandles.size());
    std::vector<Array<DeviceCluster>> compressedMeshClusters(meshHandles.size());

    tbb::parallel_for(size_t(0), clusterSourceIndices.size(), [&](size_t index) {
        size_t handleIndex = clusterSourceIndices[index];
        BufferHandle handle = meshHandles[handleIndex];
        const Mesh &mesh = container->meshes[handle.index];

        if (index != handleIndex)
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

        MeshletClustersResult &result = clusterResults[handleIndex];
        ClusterizeMesh(container, mesh, filteredIndices, result);
        compressedMeshClusters[handleIndex].Resize(result.meshlets.size());
    });

    std::vector<uint32_t> clusterOffsets;
    const size_t totalNumClusters = CalculateBufferOffsets(compressedMeshClusters, clusterOffsets);

    container->compressedClusterBytes.clear();
    container->compressedClusterBytes.resize(totalNumClusters);
    const Timer compressClustersTimer = StartTimer();
    tbb::parallel_for(size_t(0), clusterSourceIndices.size(), [&](size_t index) {
        size_t handleIndex = clusterSourceIndices[index];
        if (index != handleIndex)
        {
            return;
        }

        YBI_ASSERT(clusterOffsets[handleIndex] != INVALID_OFFSET);
        CompressClusters(clusterResults[handleIndex],
                         container->meshes[meshHandles[handleIndex].index],
                         container,
                         compressedMeshClusters[handleIndex],
                         clusterOffsets[handleIndex]);
    });

    const double compressClustersMilliseconds = StopTimer(compressClustersTimer);
    printf("CompressClusters parallel_for: %.3f ms\n", compressClustersMilliseconds);

    std::vector<uint32_t> compressedClusterByteOffsets;
    const size_t totalCompressedClusterBytes = CalculateCompressedClusterByteOffsets(
        container->compressedClusterBytes, compressedClusterByteOffsets);

    for (size_t index = 0; index < clusterSourceIndices.size(); index++)
    {
        size_t handleIndex = clusterSourceIndices[index];
        if (index != handleIndex)
        {
            continue;
        }

        YBI_ASSERT(clusterOffsets[handleIndex] != INVALID_OFFSET);
        Array<DeviceCluster> &clusters = compressedMeshClusters[handleIndex];
        for (size_t clusterIndex = 0; clusterIndex < clusters.size(); clusterIndex++)
        {
            const size_t globalClusterIndex = clusterOffsets[handleIndex] + clusterIndex;
            YBI_ASSERT(globalClusterIndex < compressedClusterByteOffsets.size());
            clusters[clusterIndex].byteOffset = compressedClusterByteOffsets[globalClusterIndex];
        }
    }

    for (size_t handleIndex = 0; handleIndex < meshHandles.size(); ++handleIndex)
    {
        const size_t sourceHandleIndex = clusterSourceIndices[handleIndex];

        YBI_ASSERT(sourceHandleIndex < clusterOffsets.size());
        YBI_ASSERT(clusterOffsets[sourceHandleIndex] != INVALID_OFFSET);

        const size_t meshIndex = meshHandles[handleIndex].index;
        Object &object = objects[meshIndex];
        object.triangleIndexOffset = INVALID_OFFSET;
        object.clusterOffset = clusterOffsets[sourceHandleIndex];
    }

    Array<BVH> &bvhs = container->GetBVHs(GeometryType::Mesh);

    const Timer buildClusterBVHsTimer = StartTimer();
    for (size_t handleIndex = 0; handleIndex < meshHandles.size(); handleIndex++)
    {
        const size_t meshIndex = meshHandles[handleIndex].index;
        Mesh &mesh = container->meshes[meshIndex];
        const size_t sourceIndex = bvhSourceIndices[handleIndex];

        if (sourceIndex != handleIndex)
        {
            const size_t sourceMeshIndex = meshHandles[sourceIndex].index;
            YBI_ASSERT(objectBVHIndices[sourceMeshIndex] >= 0);
            objectBVHIndices[meshIndex] = objectBVHIndices[sourceMeshIndex];
        }
        else
        {
            const MeshletClustersResult &result = clusterResults[handleIndex];
            YBI_ASSERT(result.meshlets.size());
            bvhs.EmplaceBack(device->BuildClusterBVH(deviceArena, container, mesh, result));
            objectBVHIndices[meshIndex] = static_cast<int>(bvhs.size() - 1);
        }
    }
    const double buildClusterBVHsMilliseconds = StopTimer(buildClusterBVHsTimer);
    printf("BuildClusterBVHs: %.3f ms\n", buildClusterBVHsMilliseconds);

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

void GeometryManager::Update(Device *device,
                             DeviceArena *deviceArena,
                             DataContainer *container,
                             Scene *rootScene,
                             DeviceScene *deviceScene)
{
    // TODO: this function is a bit overloaded maybe
    ScratchArena scratch;

    deviceScene->positions = nullptr;
    deviceScene->normals = nullptr;
    deviceScene->uvs = nullptr;
    deviceScene->uvIndices = nullptr;
    deviceScene->triangleIndices = nullptr;
    deviceScene->triangleMaterialIndices = nullptr;
    deviceScene->compressedClusters = nullptr;
    deviceScene->compressedClusterBytes = nullptr;

    YBI_ASSERT(container->defaultRectHandle.IsValid());
    YBI_ASSERT(container->defaultRectHandle.index < container->meshes.size());
    const Mesh defaultRectMesh = container->meshes[container->defaultRectHandle.index];
    for (RectLight &rectLight : container->rectLights)
    {
        if (rectLight.meshIndex == ~0u)
        {
            rectLight.meshIndex = static_cast<uint32_t>(container->meshes.size());
            container->meshes.push_back(defaultRectMesh);

            const AffineTransform lightTransform =
                rectLight.localToWorld *
                AffineTransform::Scale(rectLight.width, rectLight.height, 1.f);
            rootScene->instances.EmplaceBack(
                lightTransform,
                GeometryHandle{GeometryType::Mesh, static_cast<int>(rectLight.meshIndex)});
        }

        YBI_ASSERT(rectLight.meshIndex < container->meshes.size());
    }

    // Set up device scene
    Array<Object> objects;
    objects.Reserve(container->meshes.size());
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

#if 0
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
        YBI_ASSERT(mesh.PrimvarsAreSorted());
        size_t numPrimvars = mesh.GetPrimvars().size();

        meshPrimvarOffsets[meshIndex] = primvarOffset;
        primvarOffset += numPrimvars;
    }

    std::vector<size_t> bvhSourceIndices;
    bvhSourceIndices.reserve(container->meshes.size());
    std::unordered_map<uint64_t, size_t> bvhPairToMeshIndex;

    for (size_t meshIndex = 0; meshIndex < container->meshes.size(); meshIndex++)
    {
        Mesh &mesh = container->meshes[meshIndex];

        YBI_ASSERT(objects.size() == meshIndex);
        objects.EmplaceBack();
        Object &object = objects.Last();
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
    }
#endif

    // Object IDs reported by the intersection API are mesh indices, so keep the
    // device object table in that same index domain regardless of which build path
    // handles each mesh.
    std::vector<uint32_t> meshPrimvarOffsets(container->meshes.size());
    size_t primvarOffset = 0;
    objects.Resize(container->meshes.size());
    for (size_t meshIndex = 0; meshIndex < container->meshes.size(); meshIndex++)
    {
        const Mesh &mesh = container->meshes[meshIndex];
        YBI_ASSERT(mesh.PrimvarsAreSorted());

        meshPrimvarOffsets[meshIndex] = static_cast<uint32_t>(primvarOffset);
        primvarOffset += mesh.GetPrimvars().size();

        Object &object = objects[meshIndex];
        object = {};
        object.positionOffset = INVALID_OFFSET;
        object.normalOffset = INVALID_OFFSET;
        object.uvOffset = INVALID_OFFSET;
        object.uvIndexOffset = INVALID_OFFSET;
        object.triangleIndexOffset = INVALID_OFFSET;
        object.triangleMaterialOffset = INVALID_OFFSET;
        object.clusterOffset = INVALID_OFFSET;
        object.uvType = mesh.GetUVType();
        object.normalType = mesh.GetNormalType();
        object.primvarOffset = meshPrimvarOffsets[meshIndex];
        object.primvarCount = static_cast<uint32_t>(mesh.GetPrimvars().size());

        const uint32_t materialIndex = mesh.GetMaterialIndex();
        if (materialIndex < container->materialHandles.size())
        {
            object.materialIndex = static_cast<int>(materialIndex);
        }
    }

    for (size_t lightIndex = 0; lightIndex < container->rectLights.size(); lightIndex++)
    {
        const uint32_t meshIndex = container->rectLights[lightIndex].meshIndex;
        YBI_ASSERT(meshIndex < objects.size());
        objects[meshIndex].lightIndex = static_cast<int>(lightIndex);
    }

    // Create BVHs
    container->GetBVHs(GeometryType::Mesh).Reserve(container->meshes.size());
    printf("num meshes: %zu\n", container->meshes.size());

    Array<BufferHandle> clusterMeshHandles;
    Array<BufferHandle> triangleMeshHandles;
    // Maps object to bvh index
    Array<int> objectToBVHIndex;
    objectToBVHIndex.Resize(container->meshes.size());
    for (size_t meshIndex = 0; meshIndex < objectToBVHIndex.size(); meshIndex++)
    {
        objectToBVHIndex[meshIndex] = -1;
    }
#ifdef YBI_USE_CLUSTERS
    if (device->SupportsClusterBVHs())
    {
        for (uint32_t i = 0; i < container->meshes.size(); i++)
        {
            if (i == container->defaultRectHandle.index || objects[i].lightIndex >= 0)
            {
                continue;
            }

            clusterMeshHandles.EmplaceBack(BufferHandle{i});
            // TODO: more criteria for using clusters vs triangles
            // if (i == container->defaultRectHandle.index)
            // {
            //     triangleMeshHandles.EmplaceBack(BufferHandle{i});
            // }
            // else
            // {
            //     clusterMeshHandles.EmplaceBack(BufferHandle{i});
            // }
        }
    }
    else
#endif
    {
        for (uint32_t i = 0; i < container->meshes.size(); i++)
        {
            if (i == container->defaultRectHandle.index || objects[i].lightIndex >= 0)
            {
                continue;
            }

            triangleMeshHandles.EmplaceBack(BufferHandle{i});
        }
    }

    if (clusterMeshHandles.size())
    {
        BuildClusterMeshBVHsAndUpload(device,
                                      deviceArena,
                                      container,
                                      deviceScene,
                                      clusterMeshHandles,
                                      objects,
                                      objectToBVHIndex);
    }
    if (triangleMeshHandles.size())
    {
        BuildTriangleMeshBVHsAndUpload(device,
                                       deviceArena,
                                       container,
                                       deviceScene,
                                       triangleMeshHandles,
                                       objects,
                                       objectToBVHIndex);
    }

    if (container->rectLights.size())
    {
        Array<BVH> &meshBVHs = container->GetBVHs(GeometryType::Mesh);
        meshBVHs.EmplaceBack(container->defaultRectBVH);
        const int rectLightBVHIndex = static_cast<int>(meshBVHs.size() - 1);

        objectToBVHIndex[container->defaultRectHandle.index] = rectLightBVHIndex;
        for (const RectLight &rectLight : container->rectLights)
        {
            YBI_ASSERT(rectLight.meshIndex < objectToBVHIndex.size());
            objectToBVHIndex[rectLight.meshIndex] = rectLightBVHIndex;
        }
    }

    BuildTLAS(device, deviceArena, container, rootScene, objectToBVHIndex);

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

    // Lights
    deviceScene->deviceLights = nullptr;
    deviceScene->numLights = static_cast<uint32_t>(container->rectLights.size());
    if (deviceScene->numLights)
    {
        View<DeviceLight> hostLights =
            scratch.arena->PushArray<DeviceLight>(deviceScene->numLights);
        for (size_t lightIndex = 0; lightIndex < container->rectLights.size(); ++lightIndex)
        {
            const RectLight &rectLight = container->rectLights[lightIndex];
            DeviceLight &deviceLight = hostLights[lightIndex];

            deviceLight.worldFromLocal = rectLight.localToWorld;
            const bool invertible = rectLight.localToWorld.TryInverse(&deviceLight.localFromWorld);
            YBI_ASSERT(invertible);

            deviceLight.width = rectLight.width;
            deviceLight.height = rectLight.height;
            deviceLight.radiance = rectLight.radiance;
        }

        deviceScene->deviceLights = reinterpret_cast<DeviceLight *>(
            device->Alloc(sizeof(DeviceLight) * hostLights.size()).Data());
        device->CopyToDevice(deviceScene->deviceLights, hostLights.data(), hostLights.numBytes());
    }

    {
    }

    size_t objectBytes = sizeof(Object) * objects.size();
    deviceScene->objects = reinterpret_cast<Object *>(device->Alloc(objectBytes).Data());
    device->CopyToDevice(deviceScene->objects, objects.data(), objectBytes);

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
                                                               MaterialNodeId member,
                                                               std::vector<uint32_t> &nodes)
{
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

        mesh.SortPrimvars();
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
            ProcessMaterialGraph(material, material.base_color, shaderNodes);
        deviceMaterial.base_weight =
            ProcessMaterialGraph(material, material.base_weight, shaderNodes);
        deviceMaterial.specular_weight =
            ProcessMaterialGraph(material, material.specular_weight, shaderNodes);
        deviceMaterial.specular_color =
            ProcessMaterialGraph(material, material.specular_color, shaderNodes);
        deviceMaterial.specular_roughness =
            ProcessMaterialGraph(material, material.specular_roughness, shaderNodes);
        deviceMaterial.specular_roughness_anisotropy =
            ProcessMaterialGraph(material, material.specular_roughness_anisotropy, shaderNodes);
        deviceMaterial.specular_ior =
            ProcessMaterialGraph(material, material.specular_ior, shaderNodes);
        deviceMaterial.fuzz_weight =
            ProcessMaterialGraph(material, material.fuzz_weight, shaderNodes);
        deviceMaterial.fuzz_color =
            ProcessMaterialGraph(material, material.fuzz_color, shaderNodes);
        deviceMaterial.fuzz_roughness =
            ProcessMaterialGraph(material, material.fuzz_roughness, shaderNodes);
        // #define VISIT_PARAMETER(name, type) \
        //     deviceMaterial.name = ProcessMaterialGraph(material, #name,
        //     material.name, shaderNodes);
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
