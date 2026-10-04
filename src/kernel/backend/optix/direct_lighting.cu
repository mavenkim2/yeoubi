#pragma once

#include "kernel/backend/optix/intersection_programs.h"
#include "kernel/integrator/direct_lighting.h"
#include "kernel/launch/wavefront_params.h"
#include <optix_device.h>

namespace ybi
{

extern "C"
{
    __constant__ WavefrontParams params;
}

extern "C" __global__ void __raygen__direct_lighting()
{
    const uint32_t threadIndex = optixGetLaunchIndex().x;
    if (threadIndex >= params.numSampleIndices)
    {
        return;
    }

    uint32_t sampleIndex = params.sampleIndices[threadIndex];
    DirectLighting(params, &params.sampleState, sampleIndex);
}

} // namespace ybi
