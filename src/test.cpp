#include "device/cuda_device.h"
#include "io/usd/load.h"
#include "scene/scene.h"
#include "util/float3.h"
#include <cassert>
#include <chrono>
#include <cstdio>
#include <string>
#include <vector>

using namespace ybi;

int main(int argc, char **argv)
{
    if (argc < 2)
    {
        printf("usage: %s <usd_file>\n", argv[0]);
        return 1;
    }
    std::string usdFilePath = argv[1];

    ScenePool scenePool;
    LoadUSDScene(&scenePool, usdFilePath);
    if (scenePool.scenes.empty() || scenePool.rootSceneIndex >= scenePool.scenes.size())
    {
        fprintf(stderr, "Failed to load USD scene or invalid root: %s\n", usdFilePath.c_str());
        return 1;
    }

#if defined(WITH_CUDA) && defined(WITH_OPTIX)
    ScenePool flattenedScenePool;
    std::string flattenError;
    if (!FlattenScenePoolToRootChildren(&scenePool, &flattenedScenePool, &flattenError))
    {
        fprintf(stderr, "Failed to flatten USD ScenePool: %s\n", flattenError.c_str());
        return 1;
    }
    if (flattenedScenePool.scenes.empty() ||
        flattenedScenePool.rootSceneIndex >= flattenedScenePool.scenes.size())
    {
        fprintf(stderr, "Flattened ScenePool invalid: %s\n", usdFilePath.c_str());
        return 1;
    }

    CUDADevice device;

    std::vector<SceneMeshUploadRef> meshUploadRefs;
    CollectScenePoolMeshUploadRefs(&flattenedScenePool, &meshUploadRefs);

    printf("bvh start\n");
    auto start = std::chrono::high_resolution_clock::now();
    for (const std::unique_ptr<Scene> &scenePtr : flattenedScenePool.scenes)
    {
        Scene *scene = scenePtr.get();
        assert(scene);
        device.BuildBVH(scene);
    }
    auto end = std::chrono::high_resolution_clock::now();
    printf("bvh end\n");

    std::vector<CUdeviceptr> ownedMeshBuffers;
    ownedMeshBuffers.reserve(meshUploadRefs.size() * 2);
    for (const SceneMeshUploadRef &uploadRef : meshUploadRefs)
    {
        assert(uploadRef.mesh);
        const Mesh &mesh = *uploadRef.mesh;
        if (mesh.positions.size() == 0 || mesh.indices.size() == 0)
        {
            continue;
        }

        const size_t positionsBytes = sizeof(ybi::float3) * mesh.positions.size();
        const size_t indicesBytes = sizeof(int) * mesh.indices.size();

        CUdeviceptr positionsBuffer = device.MemAllocBytes(positionsBytes);
        CUdeviceptr indicesBuffer = device.MemAllocBytes(indicesBytes);
        device.MemcpyToDevice(positionsBuffer, mesh.positions.data(), positionsBytes);
        device.MemcpyToDevice(indicesBuffer, mesh.indices.data(), indicesBytes);

        ownedMeshBuffers.push_back(positionsBuffer);
        ownedMeshBuffers.push_back(indicesBuffer);
    }

    for (CUdeviceptr buffer : ownedMeshBuffers)
    {
        device.MemFreeBytes(buffer);
    }

    const std::chrono::duration<double, std::milli> elapsed = end - start;
    printf("time elapsed: %f ms\n", elapsed.count());
    printf("bvh size: %zi\n", device.bvhTotalAllocated);
#else
    fprintf(stderr, "CUDA+OptiX build required for BVH path.\n");
    return 1;
#endif

    return 0;
}
