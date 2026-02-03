#include "device/cuda_device.h"
#include "device/device.h"
#include "io/usd/load.h"
#include <cassert>
#include <cstdio>

using namespace ybi;

int main(int argc, char **argv)
{
    // BuildBVH();
    // TODO: hardcoded

    Scene scene;
    Test(&scene);

    CUDADevice device;

    for (Mesh &mesh : scene.meshes)
    {
        BuildBVH(&device, &scene.bvh, &mesh);
    }

    printf("bvh size: %zi\n", device.bvhTotalAllocated);

    return 0;
}
