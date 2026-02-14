#include "device/cuda_device.h"
#include "device/device.h"
#include "io/usd/load.h"
#include "scene/subdivision.h"
#include "scene/subdivision_mesh.h"
#include <cassert>
#include <chrono>
#include <cstdio>

using namespace ybi;

int main(int argc, char **argv)
{
    // BuildBVH();
    // TODO: hardcoded

    // SubdivisionMesh mesh;
    // Subdivision(mesh);

    CUDADevice device;
    device.CreateGridClusterTemplates();

    Scene scene;
    Test(&scene);

    for (SubdivisionMesh &mesh : scene.subdivisionMeshes)
    {
        Subdivision(&scene, mesh);
    }


    uint32_t totalNumTriangles = 0;

    printf("bvh start\n");
    auto start = std::chrono::high_resolution_clock::now();
    BuildBVH(&device, &scene);

    // for (Mesh &mesh : scene.meshes)
    // {
    //     // BuildBVH(&device, &scene.bvh, &mesh);
    //     totalNumTriangles += mesh.numIndices / 3;
    // }
    auto end = std::chrono::high_resolution_clock::now();
    printf("bvh end\n");
    std::chrono::duration<double, std::milli> elapsed = end - start;

    printf("time elapsed: %f\n", elapsed.count());
    printf("num tris: %i\n", totalNumTriangles);
    printf("bvh size: %zi\n", device.bvhTotalAllocated);

    return 0;
}
