#pragma once

#include "kernel/backend/cuda/atomic.h"
#include "kernel/backend/cuda/warp_reduce.h"
#include "kernel/guiding/types.h"
#include <cstdio>

#if defined(__CUDACC__)

namespace ybi
{

struct CudaVMMBackend
{
    using Float = float;
    using Vector3 = ybi::Vector3<Float>;
    using SampleData = VMMSampleData<Float>;
    using VMMStatistics = ybi::VMMStatistics<Float>;
    static constexpr uint32_t SampleWidth = 1;

    YBI_GPU static uint32_t GetThreadIndex()
    {
        return threadIdx.x;
    }

    YBI_GPU static uint32_t GetThreadStride()
    {
        return blockDim.x;
    }

    YBI_GPU static void Sync()
    {
        __syncthreads();
    }

    YBI_GPU static ExclusiveScanResult ExclusiveScanAdd(uint32_t value)
    {
        __shared__ uint32_t warpOffsets[32];
        __shared__ uint32_t total;
        const uint32_t laneIndex = threadIdx.x & 31u;
        const uint32_t warpIndex = threadIdx.x >> 5;
        const uint32_t numWarps = blockDim.x >> 5;
        uint32_t inclusive = value;
#pragma unroll
        for (uint32_t offset = 1; offset < 32; offset <<= 1)
        {
            const uint32_t other = __shfl_up_sync(0xffffffffu, inclusive, offset);
            if (laneIndex >= offset)
                inclusive += other;
        }
        if (laneIndex == 31)
            warpOffsets[warpIndex] = inclusive;
        __syncthreads();

        if (warpIndex == 0)
        {
            const uint32_t warpSum = laneIndex < numWarps ? warpOffsets[laneIndex] : 0;
            uint32_t warpInclusive = warpSum;
#pragma unroll
            for (uint32_t offset = 1; offset < 32; offset <<= 1)
            {
                const uint32_t other = __shfl_up_sync(0xffffffffu, warpInclusive, offset);
                if (laneIndex >= offset)
                    warpInclusive += other;
            }
            if (laneIndex < numWarps)
                warpOffsets[laneIndex] = warpInclusive - warpSum;
            if (laneIndex == numWarps - 1)
                total = warpInclusive;
        }
        __syncthreads();
        const ExclusiveScanResult result{warpOffsets[warpIndex] + inclusive - value, total};
        __syncthreads();
        return result;
    }

    YBI_GPU static void ReduceAdd(float value, float &destination)
    {
        const float sum = WarpReduceSum(value);
        if ((threadIdx.x & 31u) == 0)
        {
            atomicAdd(&destination, sum);
        }
    }

    YBI_GPU static void ReduceMin(float value, float &destination)
    {
        const float minimum = WarpReduceMin(value);
        if ((threadIdx.x & 31u) == 0)
            AtomicMin(&destination, minimum);
    }

    YBI_GPU static void ReduceMax(float value, float &destination)
    {
        const float maximum = WarpReduceMax(value);
        if ((threadIdx.x & 31u) == 0)
            AtomicMax(&destination, maximum);
    }

    YBI_GPU static void ReduceMin(const Vec3f &value, Vec3f &destination)
    {
        ReduceMin(value.x, destination.x);
        ReduceMin(value.y, destination.y);
        ReduceMin(value.z, destination.z);
    }

    YBI_GPU static void ReduceMax(const Vec3f &value, Vec3f &destination)
    {
        ReduceMax(value.x, destination.x);
        ReduceMax(value.y, destination.y);
        ReduceMax(value.z, destination.z);
    }

    YBI_GPU static void ReduceAdd(const Vec3f &value, Vec3f &destination)
    {
        ReduceAdd(value.x, destination.x);
        ReduceAdd(value.y, destination.y);
        ReduceAdd(value.z, destination.z);
    }

    YBI_GPU static void ReduceAdd(int64_t value, int64_t &destination)
    {
        const int64_t sum = WarpReduceSum(value);
        if ((threadIdx.x & 31u) == 0)
        {
            // CUDA exposes unsigned 64-bit atomicAdd. Modular addition also
            // accumulates negative signed sums correctly when the result fits.
            atomicAdd(reinterpret_cast<unsigned long long *>(&destination),
                      static_cast<unsigned long long>(sum));
        }
    }

    YBI_GPU static void ReduceAdd(uint32_t value, uint32_t &destination)
    {
        const uint32_t sum = static_cast<uint32_t>(WarpReduceSum(static_cast<int64_t>(value)));
        if ((threadIdx.x & 31u) == 0)
            atomicAdd(&destination, sum);
    }

    YBI_GPU static void ReportIteration(uint32_t iteration,
                                        float logLikelihood,
                                        float previousLogLikelihood,
                                        bool enabled)
    {
        if (enabled && GetThreadIndex() == 0)
        {
            if (iteration == 0)
            {
                printf("VMM %u | iteration %03u | log likelihood % .6f | delta: --\n",
                       blockIdx.x,
                       iteration + 1,
                       logLikelihood);
            }
            else
            {
                printf("VMM %u | iteration %03u | log likelihood % .6f | delta %+.6f\n",
                       blockIdx.x,
                       iteration + 1,
                       logLikelihood,
                       logLikelihood - previousLogLikelihood);
            }
        }
    }
};

} // namespace ybi

#endif
