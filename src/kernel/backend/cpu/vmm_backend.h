#pragma once

#include "core/simd_float4.h"
#include "kernel/guiding/guiding.h"

#if defined(SIMD_4)
#include <cstdio>

namespace ybi
{

// One CPU worker owns each VMM. SIMD lanes process independent samples.
struct CpuSIMDVMMBackend
{
    using Float = Float4;
    using Vector3 = ybi::Vector3<Float>;
    using SampleData = VMMSampleData<Float>;
    using VMMStatistics = ybi::VMMStatistics<Float>;
    static constexpr uint32_t SampleWidth = 4;
    static uint32_t GetThreadIndex()
    {
        return 0;
    }

    static uint32_t GetThreadStride()
    {
        return 1;
    }

    static void Sync() {}

    static void ReduceAdd(Float value, float &destination)
    {
        destination += ReduceSum(value);
    }

    static void ReduceAdd(const Vector3 &value, Vec3f &destination)
    {
        ReduceAdd(value.x, destination.x);
        ReduceAdd(value.y, destination.y);
        ReduceAdd(value.z, destination.z);
    }

    static void ReportIteration(uint32_t iteration,
                                float logLikelihood,
                                float previousLogLikelihood,
                                bool enabled)
    {
        if (!enabled)
        {
            return;
        }
        if (iteration == 0)
        {
            std::printf("iteration %03u | log likelihood % .6f | delta: --\n",
                        iteration + 1,
                        logLikelihood);
        }
        else
        {
            std::printf("iteration %03u | log likelihood % .6f | delta %+.6f\n",
                        iteration + 1,
                        logLikelihood,
                        logLikelihood - previousLogLikelihood);
        }
    }
};

} // namespace ybi

#endif
