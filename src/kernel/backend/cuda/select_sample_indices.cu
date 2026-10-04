#include "kernel/state/types.h"

namespace ybi
{

// The output buffer must hold numSampleStates indices; the caller must zero the counter.
extern "C" __global__ void SelectSampleIndices(const KernelType *kernelTypes,
                                               uint32_t numSampleStates,
                                               KernelType kernelType,
                                               uint32_t *sampleIndices,
                                               uint32_t *numSampleIndices)
{
    const uint32_t sampleIndex = blockDim.x * blockIdx.x + threadIdx.x;
    if (sampleIndex >= numSampleStates || kernelTypes[sampleIndex] != kernelType)
        return;

    const uint32_t outputIndex = atomicAdd(numSampleIndices, 1u);
    sampleIndices[outputIndex] = sampleIndex;
}

} // namespace ybi
