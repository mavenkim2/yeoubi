#pragma once

#include "kernel/backend/cuda/warp_reduce.h"
#include "kernel/guiding/guiding.h"
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

    // All threads call this in the same order. Clear destinations before use,
    // then call Sync() after the reductions and before reading the results.
    YBI_GPU static void ReduceAdd(float value, float &destination)
    {
        const float sum = WarpReduceSum(value);
        if ((threadIdx.x & 31u) == 0)
        {
            atomicAdd(&destination, sum);
        }
    }

    YBI_GPU static void ReduceAdd(const Vec3f &value, Vec3f &destination)
    {
        ReduceAdd(value.x, destination.x);
        ReduceAdd(value.y, destination.y);
        ReduceAdd(value.z, destination.z);
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
