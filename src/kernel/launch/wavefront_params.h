#pragma once

#include "core/vec2u.h"
#include "kernel/launch/kernel_params.h"
#include "kernel/state/types.h"

namespace ybi
{

struct WavefrontParams : KernelParams
{
    WavefrontSampleState sampleState;
    uint32_t *numSampleStates;
    uint32_t maxNumSampleState;

    uint32_t *sampleIndices;
    uint32_t numSampleIndices;

    Vec2u outputSize;
    uint64_t cameraWorkOffset;
    uint32_t cameraSampleStateOffset;
};

} // namespace ybi
