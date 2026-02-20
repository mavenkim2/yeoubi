#include "device/cuda_device.h"
#include "device/device.h"
#include "io/usd/load.h"
#include "scene/clusterizer.h"
#include "scene/subdivision.h"
#include "scene/subdivision_mesh.h"
#include <cassert>
#include <chrono>
#include <cstdio>
#include <string>

using namespace ybi;

int main(int argc, char **argv)
{
    if (argc < 2)
    {
        printf("usage: %s <usd_file>\n", argv[0]);
        return 1;
    }
    std::string usdFilePath = argv[1];

    // BuildBVH();

    // SubdivisionMesh mesh;
    // Subdivision(mesh);

    // CUDADevice device;
    // device.CreateGridClusterTemplates();

    ScenePool scenePool;
    LoadUSDScene(&scenePool, usdFilePath);

#if 0
    // for (SubdivisionMesh &mesh : scene.subdivisionMeshes)
    // {
    //     Subdivision(&scene, mesh);
    // }

    uint32_t totalNumTriangles = 0;
    Scene &scene = *scenePool.scenes[scenePool.rootSceneIndex];

    printf("bvh start\n");
    auto start = std::chrono::high_resolution_clock::now();
    device.BuildBVH(&scenePool);

    for (Mesh &mesh : scene.meshes)
    {
        // BuildBVH(&device, &scene.bvh, &mesh);
        totalNumTriangles += mesh.indices.size() / 3;
    }
    auto end = std::chrono::high_resolution_clock::now();
    printf("bvh end\n");
    std::chrono::duration<double, std::milli> elapsed = end - start;

    printf("time elapsed: %f\n", elapsed.count());
    printf("num tris: %i\n", totalNumTriangles);
    printf("bvh size: %zi\n", device.bvhTotalAllocated);
#endif

    return 0;
}
