#include "kernel/state/types.h"

namespace ybi
{

// The caller must zero the device counts before launching and copy them back afterward.
extern "C" __global__ void CountKernels(const KernelType *kernelTypes,
                                        uint32_t numSampleStates,
                                        KernelCounts *kernelCounts)
{
    const uint32_t sampleIndex = blockDim.x * blockIdx.x + threadIdx.x;
    if (sampleIndex >= numSampleStates)
        return;

    atomicAdd(&kernelCounts->counts[kernelTypes[sampleIndex]], 1u);
}

} // namespace ybi
