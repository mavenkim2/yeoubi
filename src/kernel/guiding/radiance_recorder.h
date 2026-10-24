#pragma once

#include "core/assert.h"
#include "core/atomic.h"
#include "core/octahedral.h"
#include "core/vec3f.h"
#include "kernel/bsdf/sample.h"
#include "kernel/guiding/types.h"

namespace ybi
{

#define YBI_PATH_GUIDING_CAMERA_VERTEX -1
#define YBI_PATH_GUIDING_MIN_TRAINING_ROUGHNESS 0.01f
#define YBI_PATH_GUIDING_DISTANCE_PROPAGATION_ROUGHNESS 0.3f

YBI_CPU_GPU bool IsPathGuidingTrainingSampleEligible(uint32_t flags, float roughness)
{
    return !(flags & BSDF_DELTA) && roughness > YBI_PATH_GUIDING_MIN_TRAINING_ROUGHNESS;
}

struct RadianceRecord
{
    Vec3f position;
    uint32_t octahedralNormal;
    uint32_t octahedralDirection;
    int parentIndex;
    uint32_t references;
    Vec3f accumulatedEmission;
    Vec3f throughput;
    float pdf;
    float distance;
    float distanceScale;
    bool hasUnrecordedContinuation;
    bool trainingSampleEligible;

    YBI_CPU_GPU void InitializeIndirectSample(const Vec3f &samplePosition,
                                              const Vec3f &normal,
                                              const Vec3f &direction,
                                              int parentRecordIndex,
                                              const Vec3f &sampleThroughput,
                                              float samplePdf,
                                              float sampleDistanceScale = 0.f,
                                              bool sampleTrainingEligible = true)
    {
        position = samplePosition;
        octahedralNormal = EncodeOctahedral(normal);
        octahedralDirection = EncodeOctahedral(direction);
        parentIndex = parentRecordIndex;
        references = 1u;
        accumulatedEmission = Vec3f(0.f);
        throughput = sampleThroughput;
        pdf = samplePdf;
        distance = FLOAT_INFINITY;
        distanceScale = sampleDistanceScale;
        hasUnrecordedContinuation = false;
        trainingSampleEligible = sampleTrainingEligible;
    }

    YBI_CPU_GPU void InitializeDirectSample(const Vec3f &samplePosition,
                                            const Vec3f &normal,
                                            const Vec3f &direction,
                                            int parentRecordIndex,
                                            const Vec3f &sampleThroughput,
                                            float samplePdf,
                                            float sampleDistance,
                                            const Vec3f &sampleAccumulatedEmission,
                                            bool sampleTrainingEligible = true)
    {
        InitializeIndirectSample(samplePosition,
                                 normal,
                                 direction,
                                 parentRecordIndex,
                                 sampleThroughput,
                                 samplePdf,
                                 0.f,
                                 sampleTrainingEligible);
        distance = sampleDistance;
        accumulatedEmission = sampleAccumulatedEmission;
    }

    YBI_CPU_GPU void AtomicAddEmission(const Vec3f &emission)
    {
        AtomicAdd(&accumulatedEmission.x, emission.x, AtomicOrder::Relaxed);
        AtomicAdd(&accumulatedEmission.y, emission.y, AtomicOrder::Relaxed);
        AtomicAdd(&accumulatedEmission.z, emission.z, AtomicOrder::Relaxed);
    }
};

struct RadianceRecorder
{
    RadianceRecord *records;
    int *freeList;
    uint32_t *allocationCursor;
    uint32_t availableCount;

    int *retiredList;
    uint32_t *retiredListCount;
    uint32_t retiredListMax;

    PathGuidingTrainingSample *trainingSamples;
    uint32_t *trainingSampleCount;
    uint32_t trainingSamplesMax;

    int numBlocksAllocated;
    int maxBlocks;

    YBI_CPU_GPU RadianceRecord &GetRecord(int recordIndex) const
    {
        YBI_ASSERT(recordIndex >= 0);
        return records[recordIndex];
    }

    YBI_CPU_GPU void RecordEmission(int recordIndex, const Vec3f &contribution) const
    {
        if (recordIndex != YBI_PATH_GUIDING_CAMERA_VERTEX)
            GetRecord(recordIndex).AtomicAddEmission(contribution);
    }

    YBI_CPU_GPU void RecordDistance(int recordIndex, const Vec3f &hitPosition) const
    {
        if (recordIndex == YBI_PATH_GUIDING_CAMERA_VERTEX)
            return;

        RadianceRecord &record = GetRecord(recordIndex);
        if (record.distance != FLOAT_INFINITY || record.hasUnrecordedContinuation)
            return;

        const float distance = Length(hitPosition - record.position);
        if (distance > 0.f && IsFinite(distance))
            record.distance = distance;
    }

    YBI_CPU_GPU bool AllocateRecord(int &recordIndex) const
    {
        uint32_t freeListIndex = AtomicAdd(allocationCursor, 1u);
        if (freeListIndex >= availableCount)
        {
            return false;
        }
        recordIndex = freeList[freeListIndex];
        return true;
    }

    YBI_CPU_GPU void Free(uint32_t recordIndex) const
    {
        // Pushes to a separate list that will be combined with the free list after the current
        // kernel execution ends.
        uint32_t retiredListIndex = AtomicAdd(retiredListCount, 1u);

        if (retiredListIndex >= retiredListMax)
        {
            // TODO: ensure retired list is as long as free list
            return;
        }
        retiredList[retiredListIndex] = recordIndex;
    }

    YBI_CPU_GPU void IncrementCount(int recordIndex) const
    {
        if (recordIndex == YBI_PATH_GUIDING_CAMERA_VERTEX)
            return;
        RadianceRecord &record = GetRecord(recordIndex);
        AtomicAdd(&record.references, 1u);
    }

    YBI_CPU_GPU void MarkDone(int recordIndex) const
    {
        while (recordIndex != YBI_PATH_GUIDING_CAMERA_VERTEX)
        {
            RadianceRecord &record = GetRecord(recordIndex);
            if (AtomicSub(&record.references, 1u, AtomicOrder::AcqRel) != 1u)
                return;

            const Vec3f incidentRadiance(
                SafeDivide(record.accumulatedEmission.x, record.throughput.x),
                SafeDivide(record.accumulatedEmission.y, record.throughput.y),
                SafeDivide(record.accumulatedEmission.z, record.throughput.z));
            WritePathGuidingTrainingSample(record, incidentRadiance.MaxComponent());
            const int parentIndex = record.parentIndex;
            RecordEmission(parentIndex, record.accumulatedEmission);
            if (parentIndex != YBI_PATH_GUIDING_CAMERA_VERTEX && record.distanceScale > 0.f)
            {
                RadianceRecord &parent = GetRecord(parentIndex);
                if (!parent.hasUnrecordedContinuation)
                    parent.distance += record.distanceScale * record.distance;
            }

            Free(static_cast<uint32_t>(recordIndex));
            recordIndex = parentIndex;
        }
    }

    YBI_CPU_GPU void RecordDirectSample(const Vec3f &position,
                                        const Vec3f &normal,
                                        const Vec3f &direction,
                                        int parentRecordIndex,
                                        const Vec3f &throughput,
                                        float pdf,
                                        float distance,
                                        const Vec3f &contribution,
                                        const BsdfEvaluationProperties &properties) const
    {
        int childRecordIndex;
        if (!AllocateRecord(childRecordIndex))
        {
            RecordEmission(parentRecordIndex, contribution);
            return;
        }

        GetRecord(childRecordIndex)
            .InitializeDirectSample(position,
                                    normal,
                                    direction,
                                    parentRecordIndex,
                                    throughput,
                                    pdf,
                                    distance,
                                    contribution,
                                    IsPathGuidingTrainingSampleEligible(properties.flags,
                                                                       properties.maxRoughness));
        IncrementCount(parentRecordIndex);
        MarkDone(childRecordIndex);
    }

    YBI_CPU_GPU int RecordIndirectSample(const Vec3f &position,
                                         const Vec3f &normal,
                                         const BsdfSample &sample,
                                         int parentRecordIndex,
                                         const Vec3f &throughput) const
    {
        int childRecordIndex;
        if (!AllocateRecord(childRecordIndex))
        {
            // TODO: Revisit distance extension across skipped records.
            if (parentRecordIndex != YBI_PATH_GUIDING_CAMERA_VERTEX)
                GetRecord(parentRecordIndex).hasUnrecordedContinuation = true;
            return parentRecordIndex;
        }

        const bool propagateDistance =
            !(sample.flags & BSDF_DIFFUSE) &&
            ((sample.flags & BSDF_DELTA) ||
             ((sample.flags & BSDF_GLOSSY) && sample.sampledRoughness >= 0.f &&
              sample.sampledRoughness < YBI_PATH_GUIDING_DISTANCE_PROPAGATION_ROUGHNESS));
        const float distanceScale = propagateDistance ? sample.refractionDistanceScale : 0.f;

        GetRecord(childRecordIndex)
            .InitializeIndirectSample(position,
                                      normal,
                                      sample.direction,
                                      parentRecordIndex,
                                      throughput,
                                      sample.stochasticPdf,
                                      distanceScale,
                                      IsPathGuidingTrainingSampleEligible(sample.flags,
                                                                         sample.sampledRoughness));
        IncrementCount(parentRecordIndex);
        MarkDone(parentRecordIndex);
        return childRecordIndex;
    }

    YBI_CPU_GPU bool WritePathGuidingTrainingSample(const RadianceRecord &record,
                                                    float incidentRadiance) const
    {
        if (!record.trainingSampleEligible)
            return false;

        const uint32_t sampleIndex = AtomicAdd(trainingSampleCount, 1u);
        if (sampleIndex >= trainingSamplesMax)
            return false;

        PathGuidingTrainingSample &sample = trainingSamples[sampleIndex];
        sample.position = record.position;
        sample.octahedralDirection = record.octahedralDirection;
        sample.incidentRadiance = incidentRadiance;
        sample.pdf = record.pdf;
        sample.distance = record.distance;
        return true;
    }
};

template <typename Backend>
YBI_GPU void ReturnRetiredRecords(RadianceRecorder &recorder)
{
    const uint32_t consumedCount = Min(*recorder.allocationCursor, recorder.availableCount);
    const uint32_t retiredCount = Min(*recorder.retiredListCount, recorder.retiredListMax);
    YBI_ASSERT(*recorder.retiredListCount <= recorder.retiredListMax);
    YBI_ASSERT(retiredCount <= consumedCount);
    const uint32_t firstFreeIndex = consumedCount - retiredCount;

    for (uint32_t i = Backend::GetThreadIndex(); i < retiredCount; i += Backend::GetThreadStride())
    {
        recorder.freeList[firstFreeIndex + i] = recorder.retiredList[i];
    }
    Backend::Sync();

    if (Backend::GetThreadIndex() == 0)
    {
        *recorder.allocationCursor = firstFreeIndex;
        *recorder.retiredListCount = 0;
    }
}

} // namespace ybi
