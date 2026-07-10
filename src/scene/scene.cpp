#include "scene/scene.h"
#include "core/arena.h"
#include "core/octahedral.h"
#include "core/thread_context.h"
#include "core/vec2f.h"
#include "device/device.h"
#include <cstdint>
#include <cstring>
#include <memory>
#include <tbb/parallel_for.h>
#include <unordered_map>
#include <vector>

namespace ybi
{

Instance::Instance(AffineTransform transform, GeometryHandle handle)
    : transform(transform), handle(handle)
{
}

Mesh::Mesh(BufferHandle positions,
           BufferHandle normals,
           AttributeType normalType,
           BufferHandle uvs,
           BufferHandle uvIndices,
           AttributeType uvType,
           BufferHandle indices)
    : positions(positions), normals(normals), uvs(uvs), uvIndices(uvIndices), indices(indices),
      uvType(uvType), normalType(normalType)
{
}

BufferHandle Mesh::GetPositionsHandle() const
{
    return positions;
}

BufferHandle Mesh::GetNormalsHandle() const
{
    return normals;
}

BufferHandle Mesh::GetUVsHandle() const
{
    return uvs;
}

BufferHandle Mesh::GetUVIndicesHandle() const
{
    return uvIndices;
}

BufferHandle Mesh::GetIndicesHandle() const
{
    return indices;
}

BufferHandle Mesh::GetFaceMaterialIndicesHandle() const
{
    return faceMaterialIndices;
}

AttributeType Mesh::GetUVType() const
{
    return uvType;
}

AttributeType Mesh::GetNormalType() const
{
    return normalType;
}

uint32_t Mesh::GetMaterialIndex() const
{
    return materialIndex;
}

void Mesh::SetPositionsHandle(BufferHandle handle)
{
    positions = handle;
}

void Mesh::SetNormalsHandle(BufferHandle handle)
{
    normals = handle;
}

void Mesh::SetUVsHandle(BufferHandle handle)
{
    uvs = handle;
}

void Mesh::SetUVIndicesHandle(BufferHandle handle)
{
    uvIndices = handle;
}

void Mesh::SetIndicesHandle(BufferHandle handle)
{
    indices = handle;
}

void Mesh::SetFaceMaterialIndicesHandle(BufferHandle handle)
{
    faceMaterialIndices = handle;
}

void Mesh::SetMaterialIndex(uint32_t index)
{
    materialIndex = index;
}

// Data container
Array<BVH> &DataContainer::GetBVHs(GeometryType type)
{
    return bvhHandles[(int)type];
}

// Scene
Scene::Scene(DataContainer *container) : container(container) {}

RenderSystem::RenderSystem(Device *device, DataContainer *container, Scene *rootScene)
    : device(device), dataContainer(container), rootScene(rootScene)
{
    deviceArena = device->CreateArena();
}

void RenderSystem::Update()
{
    YBI_ASSERT(dataContainer);
    YBI_ASSERT(rootScene);
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

    devicePtr = reinterpret_cast<T *>(device->Alloc(sizeof(T) * totalSize).Data());
    for (size_t index = 0; index < buffers.size(); index++)
    {
        const Array<T> &buffer = buffers[index];
        if (buffer.size())
        {
            device->CopyToDevice(devicePtr + offsets[index],
                                 const_cast<T *>(buffer.data()),
                                 sizeof(T) * buffer.size());
        }
    }
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
    if (container->materials.size())
    {
        hostMaterials = scratch.arena->PushArray<DeviceMaterial>(container->materials.size());
    }

    // Prepare device materials
    for (size_t materialIndex = 0; materialIndex < container->materials.size(); materialIndex++)
    {
        const Material &material = container->materials[materialIndex];
        DeviceMaterial &deviceMaterial = hostMaterials[materialIndex];
        deviceMaterial.udimBase = material.udimBase;
        deviceMaterial.udimTableOffset = material.udimTableOffset;
    }

    std::vector<uint32_t> positionOffsets;
    std::vector<uint32_t> normalOffsets;
    std::vector<uint32_t> uvOffsets;
    std::vector<uint32_t> uvIndexOffsets;
    std::vector<uint32_t> triangleIndexOffsets;
    std::vector<uint32_t> triangleMaterialIndexOffsets;

    size_t totalNumVertices = CalculateBufferOffsets(container->meshPositions, positionOffsets);
    size_t totalNumNormals = CalculateBufferOffsets(container->meshNormals, normalOffsets);
    size_t totalNumUVs = CalculateBufferOffsets(container->meshUVs, uvOffsets);
    size_t totalNumUVIndices = CalculateBufferOffsets(container->meshUVIndices, uvIndexOffsets);
    size_t totalNumTriangleIndices =
        CalculateBufferOffsets(container->meshTriangleIndices, triangleIndexOffsets);
    size_t totalNumTriangleMaterialIndices =
        CalculateBufferOffsets(container->meshFaceMaterialIndices, triangleMaterialIndexOffsets);

    std::vector<size_t> bvhSourceIndices;
    bvhSourceIndices.reserve(container->meshes.size());
    std::unordered_map<uint64_t, size_t> bvhPairToMeshIndex;

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
        object.uvType = mesh.GetUVType();
        object.normalType = mesh.GetNormalType();
        object.materialIndex = -1;

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
        if (materialIndex < container->materials.size())
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

    UploadDeviceBuffers(device,
                        deviceScene->positions,
                        container->meshPositions,
                        positionOffsets,
                        totalNumVertices);
    UploadEncodedNormalBuffers(
        device, deviceScene->normals, container->meshNormals, normalOffsets, totalNumNormals);
    UploadDeviceBuffers(device, deviceScene->float2s, container->meshUVs, uvOffsets, totalNumUVs);
    UploadDeviceBuffers(
        device, deviceScene->uvIndices, container->meshUVIndices, uvIndexOffsets, totalNumUVIndices);
    UploadDeviceBuffers(device,
                        deviceScene->triangleIndices,
                        container->meshTriangleIndices,
                        triangleIndexOffsets,
                        totalNumTriangleIndices);

    deviceScene->materials = nullptr;
    if (hostMaterials.size())
    {
        deviceScene->materials = reinterpret_cast<DeviceMaterial *>(
            device->Alloc(sizeof(DeviceMaterial) * hostMaterials.size()).Data());
        device->CopyToDevice(
            deviceScene->materials, hostMaterials.data(), hostMaterials.numBytes());
    }

    UploadDeviceBuffers(device,
                        deviceScene->triangleMaterialIndices,
                        container->meshFaceMaterialIndices,
                        triangleMaterialIndexOffsets,
                        totalNumTriangleMaterialIndices);

    deviceScene->objects = reinterpret_cast<Object *>(
        device->Alloc(sizeof(Object) * container->meshes.size()).Data());
    device->CopyToDevice(deviceScene->objects, objects.data(), objects.numBytes());

    // Create BVHs
    container->GetBVHs(GeometryType::Mesh).Reserve(container->meshes.size());
    printf("num meshes: %zu\n", container->meshes.size());

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

    printf("total num vertices: %zu\n", totalNumVertices);
    printf("num instances: %zu\n", rootScene->instances.size());

    // TODO: supports 2+ level scenes
    container->GetBVHs(GeometryType::Scene)
        .EmplaceBack(device->BuildInstancesBVH(deviceArena, rootScene->instances, container));
}

} // namespace ybi
