#include "scene/scene.h"
#include "core/arena.h"
#include "core/hash.h"
#include "device/device.h"
#include <memory>
#include <unordered_map>

namespace ybi
{

Instance::Instance(AffineTransform transform, GeometryHandle handle)
    : transform(transform), handle(handle)
{
}

Mesh::Mesh(Array<Vec3f> &&positions) : positions(std::move(positions)) {}

Mesh::Mesh(Array<Vec3f> &&positions, Array<int> &&indices)
    : positions(std::move(positions)), indices(std::move(indices))
{
}

size_t Mesh::GetNumVertices() const
{
    return positions.size();
}

const Array<Vec3f> &Mesh::GetPositions() const
{
    return positions;
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
    geometryManager.Update(device, deviceArena.get(), dataContainer, rootScene);
}

void GeometryManager::Update(Device *device,
                             DeviceArena *deviceArena,
                             DataContainer *container,
                             Scene *rootScene)
{

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

// you build the BVHs <-- done
// set up optix programs
// set up optix kernels
// trace primary visibility and output some random color
// then set up cross with vulkan
// add a basic camera
// add basic input, probably glfw?
// interactivity
// move
// die

} // namespace ybi
