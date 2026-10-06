#include "kernel/state/types.h"
#include <cassert>

namespace ybi
{

// Launch with 256 threads in one dimension. The output buffer must hold
// numSampleStates indices; the caller must zero the counter. Each block writes
// one contiguous range in sample order; ranges across blocks are unordered.
extern "C" __global__ void SelectSampleIndices(const KernelType *kernelTypes,
                                               uint32_t numSampleStates,
                                               KernelType kernelType,
                                               uint32_t *sampleIndices,
                                               uint32_t *numSampleIndices)
{
    assert(blockDim.x == 256 && blockDim.y == 1 && blockDim.z == 1);
    constexpr uint32_t warpsPerBlock = 8;
    __shared__ uint32_t warpOffsets[warpsPerBlock];
    __shared__ uint32_t blockOffset;

    const uint32_t sampleIndex = blockDim.x * blockIdx.x + threadIdx.x;
    const bool predicate = sampleIndex < numSampleStates && kernelTypes[sampleIndex] == kernelType;
    const uint32_t laneIndex = threadIdx.x & 31u;
    const uint32_t warpIndex = threadIdx.x >> 5u;
    const unsigned matches = __ballot_sync(0xffffffffu, predicate);

    if (laneIndex == 0)
        warpOffsets[warpIndex] = __popc(matches);

    __syncthreads();

    if (threadIdx.x == 0)
    {
        uint32_t count = 0;
        for (uint32_t warp = 0; warp < warpsPerBlock; ++warp)
        {
            const uint32_t warpCount = warpOffsets[warp];
            warpOffsets[warp] = count;
            count += warpCount;
        }
        blockOffset = count > 0 ? atomicAdd(numSampleIndices, count) : 0;
    }

    __syncthreads();

    if (predicate)
    {
        const uint32_t rank = __popc(matches & ((1u << laneIndex) - 1u));
        const uint32_t outputIndex = blockOffset + warpOffsets[warpIndex] + rank;
        sampleIndices[outputIndex] = sampleIndex;
    }
}

} // namespace ybi
