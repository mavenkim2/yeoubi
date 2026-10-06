#pragma once

#include "core/simd_float4.h"
#include "kernel/guiding/guiding.h"
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
    uint32_t mixtureIndex{0};

    uint32_t GetThreadIndex() const
    {
        return 0;
    }

    uint32_t GetThreadStride() const
    {
        return 1;
    }

    void Sync() const {}

    void ReduceAdd(Float value, float &destination) const
    {
        destination += ReduceSum(value);
    }

    void ReduceAdd(const Vector3 &value, Vec3f &destination) const
    {
        ReduceAdd(value.x, destination.x);
        ReduceAdd(value.y, destination.y);
        ReduceAdd(value.z, destination.z);
    }

    void ReportIteration(uint32_t iteration,
                         float logLikelihood,
                         float previousLogLikelihood,
                         bool enabled) const
    {
        if (!enabled)
        {
            return;
        }
        if (iteration == 0)
        {
            std::printf("VMM %u | iteration %03u | log likelihood % .6f | delta: --\n",
                        mixtureIndex, iteration + 1, logLikelihood);
        }
        else
        {
            std::printf("VMM %u | iteration %03u | log likelihood % .6f | delta %+.6f\n",
                        mixtureIndex, iteration + 1, logLikelihood,
                        logLikelihood - previousLogLikelihood);
        }
    }
};

} // namespace ybi
