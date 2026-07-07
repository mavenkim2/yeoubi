#include "scene/scene.h"
#include "core/arena.h"
#include "core/hash.h"
#include "core/thread_context.h"
#include "device/device.h"
#include <memory>
#include <unordered_map>

namespace ybi
{

Instance::Instance(AffineTransform transform, GeometryHandle handle)
    : transform(transform), handle(handle)
{
}

Mesh::Mesh(Array<Vec3f> &&positions,
           Array<Vec2f> &&uvs,
           AttributeType uvType,
           Array<int> &&indices)
    : positions(std::move(positions)), uvs(std::move(uvs)), uvType(uvType),
      indices(std::move(indices))
{
}

size_t Mesh::GetNumVertices() const
{
    return positions.size();
}

size_t Mesh::GetNumUVs() const
{
    return uvs.size();
}

AttributeType Mesh::GetUVType() const
{
    return uvType;
}

uint32_t Mesh::GetMaterialIndex() const
{
    return materialIndex;
}

const Array<Vec3f> &Mesh::GetPositions() const
{
    return positions;
}

const Array<Vec2f> &Mesh::GetUVs() const
{
    return uvs;
}

const Array<int> &Mesh::GetIndices() const
{
    return indices;
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
    size_t totalNumUVs = 0;
    size_t totalNumTriangleIndices = 0;
    for (size_t meshIndex = 0; meshIndex < container->meshes.size(); meshIndex++)
    {
        Mesh &mesh = container->meshes[meshIndex];
        Object &object = objects[meshIndex];
        object.uvOffset = totalNumUVs;
        object.triangleIndexOffset = static_cast<uint32_t>(totalNumTriangleIndices);
        object.uvType = mesh.GetUVType();
        object.udimBase = ~0u;
        object.udimTableOffset = ~0u;

        uint32_t materialIndex = mesh.GetMaterialIndex();
        if (materialIndex < container->materials.size())
        {
            const Material &material = container->materials[materialIndex];
            object.udimBase = material.udimBase;
            object.udimTableOffset = material.udimTableOffset;
        }

        totalNumUVs += mesh.GetNumUVs();
        totalNumTriangleIndices += mesh.GetIndices().size();
    }

    View<Vec2f> hostUvs = scratch.arena->PushArray<Vec2f>(totalNumUVs);
    View<int> hostTriangleIndices = scratch.arena->PushArray<int>(totalNumTriangleIndices);
    totalNumUVs = 0;
    totalNumTriangleIndices = 0;
    for (size_t meshIndex = 0; meshIndex < container->meshes.size(); meshIndex++)
    {
        Mesh &mesh = container->meshes[meshIndex];
        const Array<Vec2f> &uvs = mesh.GetUVs();
        const Array<int> &indices = mesh.GetIndices();
        if (uvs.size())
        {
            memcpy(hostUvs.data() + totalNumUVs, uvs.data(), sizeof(uvs[0]) * uvs.size());
            totalNumUVs += uvs.size();
        }
        if (indices.size())
        {
            memcpy(hostTriangleIndices.data() + totalNumTriangleIndices,
                   indices.data(),
                   sizeof(indices[0]) * indices.size());
            totalNumTriangleIndices += indices.size();
        }
    }

    deviceScene->float2s =
        reinterpret_cast<Vec2f *>(device->Alloc(sizeof(Vec2f) * totalNumUVs).Data());
    device->CopyToDevice(deviceScene->float2s, hostUvs.data(), hostUvs.numBytes());

    deviceScene->triangleIndices =
        reinterpret_cast<int *>(device->Alloc(sizeof(int) * totalNumTriangleIndices).Data());
    device->CopyToDevice(
        deviceScene->triangleIndices, hostTriangleIndices.data(), hostTriangleIndices.numBytes());

    deviceScene->objects = reinterpret_cast<Object *>(
        device->Alloc(sizeof(Object) * container->meshes.size()).Data());
    device->CopyToDevice(deviceScene->objects, objects.data(), objects.numBytes());

    // Create BVHs
    container->GetBVHs(GeometryType::Mesh).Reserve(container->meshes.size());
    printf("num meshes: %zu\n", container->meshes.size());

    std::unordered_multimap<uint64_t, size_t> meshMap;

    size_t totalNumVertices = 0;
    for (size_t meshIndex = 0; meshIndex < container->meshes.size(); meshIndex++)
    {
        Mesh &mesh = container->meshes[meshIndex];
        size_t meshNumVertices = mesh.GetNumVertices();
        const Array<Vec3f> &meshPositions = mesh.GetPositions();
        const Array<int> &meshIndices = mesh.GetIndices();

        uint64_t hash =
            MurmurHash64A((uint8_t *)meshPositions.data(), sizeof(Vec3f) * meshNumVertices, 0);
        hash = MurmurHash64A(
            (uint8_t *)meshIndices.data(), sizeof(meshIndices[0]) * meshIndices.size(), hash);

        size_t duplicateIndex = ~0ull;
        auto range = meshMap.equal_range(hash);
        for (auto it = range.first; it != range.second; it++)
        {
            size_t otherMeshIndex = it->second;
            Mesh &otherMesh = container->meshes[otherMeshIndex];

            const Array<Vec3f> &otherPositions = otherMesh.GetPositions();
            const Array<int> &otherIndices = otherMesh.GetIndices();

            if (otherMesh.GetNumVertices() != meshNumVertices ||
                otherIndices.size() != meshIndices.size())
                continue;

            int result = memcmp(meshPositions.data(),
                                otherPositions.data(),
                                sizeof(meshPositions[0]) * meshNumVertices);

            if (result != 0)
                continue;

            result = memcmp(meshIndices.data(),
                            otherIndices.data(),
                            sizeof(meshIndices[0]) * meshIndices.size());

            if (result == 0)
            {
                duplicateIndex = otherMeshIndex;
                break;
            }
        }

        if (duplicateIndex != ~0ull)
        {
            Array<BVH> &bvhs = container->GetBVHs(GeometryType::Mesh);
            // TODO: this is brittle. also bvhs probably shouldn't be duplicated
            bvhs.EmplaceBack(bvhs[duplicateIndex]);
        }
        else
        {
            totalNumVertices += meshNumVertices;
            meshMap.insert({hash, meshIndex});

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
