#pragma once

#include "kernel/launch/wavefront_params.h"
#include "kernel/integrator/shade_hit.h"

namespace ybi
{

extern "C" __global__ void ShadeHitKernel(WavefrontParams params)
{
    uint32_t threadIndex = blockDim.x * blockIdx.x + threadIdx.x;
    if (threadIndex >= params.numSampleIndices)
    {
        return;
    }

    uint32_t sampleIndex = params.sampleIndices[threadIndex];
    ShadeHit(GetShadingContext(params), &params.sampleState, sampleIndex);
    return;
}

} // namespace ybi
