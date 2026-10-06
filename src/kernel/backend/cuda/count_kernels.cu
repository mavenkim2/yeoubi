#include "kernel/state/types.h"
#include <cassert>

namespace ybi
{

// Launch with 256 threads in one dimension. The caller must zero the device counts
// before launching and copy them back afterward.
extern "C" __global__ void
CountKernels(const KernelType *kernelTypes, uint32_t numSampleStates, KernelCounts *kernelCounts)
{
    assert(blockDim.x == 256 && blockDim.y == 1 && blockDim.z == 1);
    constexpr uint32_t warpsPerBlock = 8;
    __shared__ uint32_t warpCounts[warpsPerBlock][KERNEL_TYPE_COUNT];

    const uint32_t sampleIndex = blockDim.x * blockIdx.x + threadIdx.x;
    const KernelType kernelType =
        sampleIndex < numSampleStates ? kernelTypes[sampleIndex] : KERNEL_TYPE_COUNT;
    const uint32_t laneIndex = threadIdx.x & 31u;
    const uint32_t warpIndex = threadIdx.x >> 5u;

#pragma unroll
    for (uint32_t type = 0; type < KERNEL_TYPE_COUNT; ++type)
    {
        const unsigned matches = __ballot_sync(0xffffffffu, kernelType == type);
        if (laneIndex == 0)
            warpCounts[warpIndex][type] = __popc(matches);
    }

    __syncthreads();

    if (threadIdx.x < KERNEL_TYPE_COUNT)
    {
        uint32_t count = 0;
#pragma unroll
        for (uint32_t warp = 0; warp < warpsPerBlock; ++warp)
            count += warpCounts[warp][threadIdx.x];
        if (count > 0)
            atomicAdd(&kernelCounts->counts[threadIdx.x], count);
    }
}

} // namespace ybi
