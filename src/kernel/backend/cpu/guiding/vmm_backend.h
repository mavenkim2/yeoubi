#pragma once

#include "core/atomic.h"
#include "kernel/backend/cpu/reduction.h"
#include "kernel/guiding/types.h"

#include <cstdio>

namespace ybi
{

// One CPU worker owns each VMM.
struct CpuVMMBackend
{
    using Float = float;
    using Vector3 = ybi::Vector3<Float>;
    using SampleData = VMMSampleData<Float>;
    using VMMStatistics = ybi::VMMStatistics<Float>;
    static constexpr uint32_t SampleWidth = 1;
    static uint32_t GetThreadIndex()
    {
        return 0;
    }

    static uint32_t GetThreadStride()
    {
        return 1;
    }

    static void Sync() {}

    template <typename T, typename Function>
    static T ReduceBlock(const Function &function)
    {
        T values[VMM_BLOCK_SIZE] = {};
        for (uint32_t thread = 0; thread < VMM_BLOCK_SIZE; ++thread)
            values[thread] = function(thread);
        return BlockReduceSum<T, VMM_BLOCK_SIZE>(values);
    }

    static ExclusiveScanResult ExclusiveScanAdd(uint32_t value) { return {0, value}; }

    static void ReduceAdd(int64_t value, int64_t &destination) { destination += value; }
    static void ReduceAdd(uint32_t value, uint32_t &destination) { destination += value; }
    static void ReduceAdd(float value, float &destination) { destination += value; }

    static void ReduceMin(float value, float &destination) { AtomicMin(&destination, value); }
    static void ReduceMax(float value, float &destination) { AtomicMax(&destination, value); }

    static void ReduceMin(const Vec3f &value, Vec3f &destination)
    {
        ReduceMin(value.x, destination.x);
        ReduceMin(value.y, destination.y);
        ReduceMin(value.z, destination.z);
    }

    static void ReduceMax(const Vec3f &value, Vec3f &destination)
    {
        ReduceMax(value.x, destination.x);
        ReduceMax(value.y, destination.y);
        ReduceMax(value.z, destination.z);
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

