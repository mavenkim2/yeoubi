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

    printf("num meshes: %lli\n", container->meshes.size());

    std::unordered_multimap<uint64_t, size_t> meshMap;

    size_t totalNumVertices = 0;
    for (size_t meshIndex = 0; meshIndex < container->meshes.size(); meshIndex++)
    {
        Mesh &mesh = container->meshes[meshIndex];
        size_t meshNumVertices = mesh.GetNumVertices();
        const Array<Vec3f> &meshPositions = mesh.GetPositions();

        uint64_t hash =
            MurmurHash64A((uint8_t *)meshPositions.data(), sizeof(Vec3f) * meshNumVertices, 0);
        // printf("%llu\n", hash);

#if 1
        size_t duplicateIndex = ~0ull;
        auto range = meshMap.equal_range(hash);
        for (auto it = range.first; it != range.second; it++)
        {
            size_t otherMeshIndex = it->second;
            Mesh &otherMesh = container->meshes[otherMeshIndex];

            const Array<Vec3f> &otherPositions = otherMesh.GetPositions();
            if (otherMesh.GetNumVertices() != meshNumVertices)
                continue;

            bool repeats = true;
            for (size_t i = 0; i < meshNumVertices; i++)
            {
                if (meshPositions[i] != otherPositions[i])
                {
                    repeats = false;
                    break;
                }
            }

            if (repeats)
            {
                duplicateIndex = otherMeshIndex;
                break;
            }
        }

        if (duplicateIndex != ~0ull)
        {
            Array<BVH> &bvhs = container->GetBVHs(GeometryType::Mesh);
            // TODO: this probably shouldn't be duplicated like this
            bvhs.EmplaceBack(bvhs[duplicateIndex]);
        }
        else
        {
            totalNumVertices += meshNumVertices;
            meshMap.insert({hash, meshIndex});

            container->GetBVHs(GeometryType::Mesh)
                .EmplaceBack(device->BuildTriangleBVH(deviceArena, mesh));
        }
#endif

        // printf("mesh: %lli, verts: %lli\n", meshIndex, mesh.GetNumVertices());
    }

    printf("total num vertices: %lli\n", totalNumVertices);
    printf("num instances: %lli\n", rootScene->instances.size());

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
