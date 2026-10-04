#include "kernel/backend/optix/intersection_programs.h"
#include "kernel/integrator/secondary_ray.h"
#include "kernel/launch/wavefront_params.h"
#include <optix_device.h>

namespace ybi
{

extern "C"
{
    __constant__ WavefrontParams params;
}

extern "C" __global__ void __raygen__secondary()
{
    const uint32_t threadIndex = optixGetLaunchIndex().x;
    if (threadIndex >= params.numSampleIndices)
    {
        return;
    }

    uint32_t sampleIndex = params.sampleIndices[threadIndex];
    SecondaryRay(params, &params.sampleState, sampleIndex);
}

} // namespace ybi
