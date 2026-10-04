#include "kernel/integrator/generate_camera_ray.h"

namespace ybi
{

// Launch in one dimension. The caller advances workOffset by numWorkItems,
// with a total of imageWidth * imageHeight work items per sample pass.
extern "C" __global__ void GenerateCameraRays(const WavefrontParams *kernelParams,
                                              uint32_t numWorkItems)
{
    const WavefrontParams &params = *kernelParams;
    YBI_ASSERT(uint64_t(params.cameraSampleStateOffset) + numWorkItems <= params.maxNumSampleState);
    YBI_ASSERT(params.cameraWorkOffset + numWorkItems <= uint64_t(params.imageWidth) * params.imageHeight);
    const uint32_t launchIndex = blockDim.x * blockIdx.x + threadIdx.x;
    if (launchIndex >= numWorkItems)
        return;

    GenerateCameraRay(params,
                       params.cameraSampleStateOffset + launchIndex,
                       params.cameraWorkOffset + launchIndex);
}

} // namespace ybi
