#include "scene/scene.h"
#include "core/arena.h"
#include "core/thread_context.h"
#include "device/device.h"
#include <cstdint>
#include <cstring>
#include <memory>
#include <vector>

namespace ybi
{

Instance::Instance(AffineTransform transform, GeometryHandle handle)
    : transform(transform), handle(handle)
{
}

Mesh::Mesh(Array<Vec3f> &&positions,
           Array<Vec3f> &&normals,
           AttributeType normalType,
           Array<Vec2f> &&uvs,
           AttributeType uvType,
           Array<int> &&indices)
    : positions(std::move(positions)), normals(std::move(normals)), uvs(std::move(uvs)),
      indices(std::move(indices)), uvType(uvType), normalType(normalType)
{
}

size_t Mesh::GetNumVertices() const
{
    return positions.size();
}

size_t Mesh::GetNumNormals() const
{
    return normals.size();
}

size_t Mesh::GetNumUVs() const
{
    return uvs.size();
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

const Array<Vec3f> &Mesh::GetPositions() const
{
    return positions;
}

const Array<Vec3f> &Mesh::GetNormals() const
{
    return normals;
}

const Array<Vec2f> &Mesh::GetUVs() const
{
    return uvs;
}

const Array<int> &Mesh::GetIndices() const
{
    return indices;
}

const Array<int> &Mesh::GetFaceMaterialIndices() const
{
    return faceMaterialIndices;
}

void Mesh::SetMaterialIndex(uint32_t index)
{
    materialIndex = index;
}

void Mesh::SetFaceMaterialIndices(Array<int> &&indices)
{
    faceMaterialIndices = std::move(indices);
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
struct UniqueArrayUpload
{
    const Array<T> *array = nullptr;
    size_t offset = 0;
};

template <typename T>
static void CopyUniqueArrays(View<T> dst, const std::vector<UniqueArrayUpload<T>> &uploads)
{
    for (const UniqueArrayUpload<T> &upload : uploads)
    {
        const Array<T> &array = *upload.array;
        if (array.size())
        {
            memcpy(dst.data() + upload.offset, array.data(), sizeof(T) * array.size());
        }
    }
}

template <typename T>
static void UploadDeviceArray(Device *device, T *&devicePtr, View<T> host)
{
    devicePtr = nullptr;
    if (host.size())
    {
        devicePtr = reinterpret_cast<T *>(device->Alloc(sizeof(T) * host.size()).Data());
        device->CopyToDevice(devicePtr, host.data(), host.numBytes());
    }
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

    // Prepare device objects
    size_t totalNumUVs = 0;
    size_t totalNumTriangleMaterialIndices = 0;
    size_t totalNumVertices = 0;
    size_t totalNumNormals = 0;
    size_t totalNumTriangleIndices = 0;

    DataCache positionCache;
    DataCache normalCache;
    DataCache uvCache;
    DataCache indicesCache;
    DataCache triangleMaterialIndexCache;

    std::vector<UniqueArrayUpload<Vec3f>> positionUploads;
    std::vector<UniqueArrayUpload<Vec3f>> normalUploads;
    std::vector<UniqueArrayUpload<Vec2f>> uvUploads;
    std::vector<UniqueArrayUpload<int>> triangleIndexUploads;
    std::vector<UniqueArrayUpload<int>> triangleMaterialIndexUploads;
    std::vector<size_t> bvhSourceIndices;
    positionUploads.reserve(container->meshes.size());
    normalUploads.reserve(container->meshes.size());
    uvUploads.reserve(container->meshes.size());
    triangleIndexUploads.reserve(container->meshes.size());
    triangleMaterialIndexUploads.reserve(container->meshes.size());
    bvhSourceIndices.reserve(container->meshes.size());

    for (size_t meshIndex = 0; meshIndex < container->meshes.size(); meshIndex++)
    {
        Mesh &mesh = container->meshes[meshIndex];
        const Array<Vec3f> &meshPositions = mesh.GetPositions();
        const Array<Vec3f> &meshNormals = mesh.GetNormals();
        const Array<Vec2f> &meshUvs = mesh.GetUVs();
        const Array<int> &meshIndices = mesh.GetIndices();
        const Array<int> &faceMaterialIndices = mesh.GetFaceMaterialIndices();

        Object &object = objects[meshIndex];
        object.positionOffset = INVALID_OFFSET;
        object.normalOffset = INVALID_OFFSET;
        object.uvOffset = INVALID_OFFSET;
        object.triangleIndexOffset = INVALID_OFFSET;
        object.triangleMaterialOffset = INVALID_OFFSET;
        object.uvType = mesh.GetUVType();
        object.normalType = mesh.GetNormalType();
        object.materialIndex = -1;

        size_t meshPositionIndex = positionCache.GetOrAddArray(
            meshPositions, meshIndex, [&](size_t existingIndex) -> const Array<Vec3f> & {
                return container->meshes[existingIndex].GetPositions();
            });

        if (meshPositions.size())
        {
            if (meshPositionIndex == INVALID_INDEX)
            {
                object.positionOffset = static_cast<uint32_t>(totalNumVertices);
                positionUploads.push_back({&meshPositions, totalNumVertices});
                totalNumVertices += meshPositions.size();
            }
            else
            {
                object.positionOffset = objects[meshPositionIndex].positionOffset;
            }
        }

        if (meshNormals.size())
        {
            size_t meshNormalIndex = normalCache.GetOrAddArray(
                meshNormals, meshIndex, [&](size_t existingIndex) -> const Array<Vec3f> & {
                    return container->meshes[existingIndex].GetNormals();
                });

            if (meshNormalIndex == INVALID_INDEX)
            {
                object.normalOffset = static_cast<uint32_t>(totalNumNormals);
                normalUploads.push_back({&meshNormals, totalNumNormals});
                totalNumNormals += meshNormals.size();
            }
            else
            {
                object.normalOffset = objects[meshNormalIndex].normalOffset;
            }
        }

        size_t meshIndicesIndex = indicesCache.GetOrAddArray(
            meshIndices, meshIndex, [&](size_t existingIndex) -> const Array<int> & {
                return container->meshes[existingIndex].GetIndices();
            });

        if (meshIndices.size())
        {
            if (meshIndicesIndex == INVALID_INDEX)
            {
                object.triangleIndexOffset = static_cast<uint32_t>(totalNumTriangleIndices);
                triangleIndexUploads.push_back({&meshIndices, totalNumTriangleIndices});
                totalNumTriangleIndices += meshIndices.size();
            }
            else
            {
                object.triangleIndexOffset = objects[meshIndicesIndex].triangleIndexOffset;
            }
        }

        if (meshUvs.size())
        {
            size_t meshUvIndex = uvCache.GetOrAddArray(
                meshUvs, meshIndex, [&](size_t existingIndex) -> const Array<Vec2f> & {
                    return container->meshes[existingIndex].GetUVs();
                });
            if (meshUvIndex == INVALID_INDEX)
            {
                object.uvOffset = static_cast<uint32_t>(totalNumUVs);
                uvUploads.push_back({&meshUvs, totalNumUVs});
                totalNumUVs += meshUvs.size();
            }
            else
            {
                object.uvOffset = objects[meshUvIndex].uvOffset;
            }
        }

        if (faceMaterialIndices.size())
        {
            YBI_ASSERT(faceMaterialIndices.size() == meshIndices.size() / 3);
            size_t faceMaterialIndex = triangleMaterialIndexCache.GetOrAddArray(
                faceMaterialIndices, meshIndex, [&](size_t existingIndex) -> const Array<int> & {
                    return container->meshes[existingIndex].GetFaceMaterialIndices();
                });

            if (faceMaterialIndex == INVALID_INDEX)
            {
                object.triangleMaterialOffset =
                    static_cast<uint32_t>(totalNumTriangleMaterialIndices);
                triangleMaterialIndexUploads.push_back(
                    {&faceMaterialIndices, totalNumTriangleMaterialIndices});
                totalNumTriangleMaterialIndices += faceMaterialIndices.size();
            }
            else
            {
                object.triangleMaterialOffset = objects[faceMaterialIndex].triangleMaterialOffset;
            }
        }

        uint32_t materialIndex = mesh.GetMaterialIndex();
        if (materialIndex < container->materials.size())
        {
            object.materialIndex = static_cast<int>(materialIndex);
        }

        if (meshPositions.size() && meshIndices.size() && meshPositionIndex != INVALID_INDEX &&
            meshIndicesIndex != INVALID_INDEX && meshPositionIndex == meshIndicesIndex)
        {
            bvhSourceIndices.push_back(meshPositionIndex);
        }
        else
        {
            bvhSourceIndices.push_back(meshIndex);
        }
    }

    View<Vec3f> hostPositions;
    View<Vec3f> hostNormals;
    View<Vec2f> hostUvs;
    View<int> hostTriangleIndices;
    View<int> hostTriangleMaterialIndices;
    if (totalNumVertices)
    {
        hostPositions = scratch.arena->PushArray<Vec3f>(totalNumVertices);
    }
    if (totalNumNormals)
    {
        hostNormals = scratch.arena->PushArray<Vec3f>(totalNumNormals);
    }
    if (totalNumUVs)
    {
        hostUvs = scratch.arena->PushArray<Vec2f>(totalNumUVs);
    }
    if (totalNumTriangleIndices)
    {
        hostTriangleIndices = scratch.arena->PushArray<int>(totalNumTriangleIndices);
    }
    if (totalNumTriangleMaterialIndices)
    {
        hostTriangleMaterialIndices =
            scratch.arena->PushArray<int>(totalNumTriangleMaterialIndices);
    }

    CopyUniqueArrays(hostPositions, positionUploads);
    CopyUniqueArrays(hostNormals, normalUploads);
    CopyUniqueArrays(hostUvs, uvUploads);
    CopyUniqueArrays(hostTriangleIndices, triangleIndexUploads);
    CopyUniqueArrays(hostTriangleMaterialIndices, triangleMaterialIndexUploads);

    UploadDeviceArray(device, deviceScene->positions, hostPositions);
    UploadDeviceArray(device, deviceScene->normals, hostNormals);
    UploadDeviceArray(device, deviceScene->float2s, hostUvs);
    UploadDeviceArray(device, deviceScene->triangleIndices, hostTriangleIndices);

    deviceScene->materials = nullptr;
    if (hostMaterials.size())
    {
        deviceScene->materials = reinterpret_cast<DeviceMaterial *>(
            device->Alloc(sizeof(DeviceMaterial) * hostMaterials.size()).Data());
        device->CopyToDevice(
            deviceScene->materials, hostMaterials.data(), hostMaterials.numBytes());
    }

    UploadDeviceArray(device, deviceScene->triangleMaterialIndices, hostTriangleMaterialIndices);

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
                .EmplaceBack(device->BuildTriangleBVH(deviceArena, mesh));
        }
    }

    printf("total num vertices: %zu\n", totalNumVertices);
    printf("num instances: %zu\n", rootScene->instances.size());

    // TODO: supports 2+ level scenes
    container->GetBVHs(GeometryType::Scene)
        .EmplaceBack(device->BuildInstancesBVH(deviceArena, rootScene->instances, container));
}

} // namespace ybi
