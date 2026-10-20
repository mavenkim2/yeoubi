#pragma once

#include "kernel/launch/kernel_params.h"
#include "kernel/state/types.h"

namespace ybi
{

template <typename SampleState>
YBI_GPU void TerminatePath(const KernelParams &params,
                           const SampleState *sampleState,
                           uint32_t sampleIndex,
                           const Vec3f &radiance)
{
    // One sample per pixel owns this output, so no atomic accumulation is needed.
    const uint32_t pixelIndex = sampleState->ReadPixelIndex(sampleIndex);
    params.kernelOutput[pixelIndex] = radiance;
#if defined(YBI_USE_PATH_GUIDING)
    params.radianceRecorder.MarkDone(sampleState->ReadRadianceRecordIndex(sampleIndex));
#endif
    sampleState->WriteKernelType(sampleIndex, KERNEL_FINISHED);
}

template <typename SampleState>
YBI_GPU void
TerminatePath(const KernelParams &params, const SampleState *sampleState, uint32_t sampleIndex)
{
    TerminatePath(params, sampleState, sampleIndex, sampleState->ReadRadiance(sampleIndex));
}

} // namespace ybi
