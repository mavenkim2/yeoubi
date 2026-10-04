#pragma once

#include "kernel/launch/kernel_params.h"
#include "kernel/state/types.h"

namespace ybi
{

struct WavefrontParams : KernelParams
{
    // Wavefront
    WavefrontSampleState sampleState;
    uint32_t *numSampleStates;
    uint32_t maxNumSampleState;

    uint32_t *sampleIndices;
    uint32_t numSampleIndices;

    uint64_t cameraWorkOffset;
    uint32_t cameraSampleStateOffset;
};

} // namespace ybi
