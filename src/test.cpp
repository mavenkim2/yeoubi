#include "device/cuda_device.h"
#include "device/device.h"
#include "io/usd/load.h"
#include <cassert>
#include <chrono>
#include <cstdio>

using namespace ybi;

int main(int argc, char **argv)
{
    // BuildBVH();
    // TODO: hardcoded

    Scene scene;
    Test(&scene);

    CUDADevice device;

    uint32_t totalNumTriangles = 0;
    auto start = std::chrono::high_resolution_clock::now();
    BuildBVH(&device, &scene);

    // for (Mesh &mesh : scene.meshes)
    // {
    //     // BuildBVH(&device, &scene.bvh, &mesh);
    //     totalNumTriangles += mesh.numIndices / 3;
    // }
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> elapsed = end - start;

    printf("time elapsed: %f\n", elapsed.count());
    printf("num tris: %i\n", totalNumTriangles);
    printf("bvh size: %zi\n", device.bvhTotalAllocated);

    return 0;
}
