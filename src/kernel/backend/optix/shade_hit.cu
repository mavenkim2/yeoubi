#pragma once

#include "kernel/integrator/shade_hit.h"
#include "kernel/launch/wavefront_params.h"
#include <optix_device.h>

namespace ybi
{

extern "C"
{
    __constant__ WavefrontParams params;
}

extern "C" __global__ void __raygen__shade_hit()
{
    const uint32_t threadIndex = optixGetLaunchIndex().x;
    if (threadIndex >= params.numSampleIndices)
    {
        return;
    }

    uint32_t sampleIndex = params.sampleIndices[threadIndex];
    ShadeHit(GetShadingContext(params),
             static_cast<const KernelParams *>(&params),
             &params.sampleState,
             sampleIndex);
}

} // namespace ybi
