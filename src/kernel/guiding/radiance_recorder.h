#pragma once

#include "core/assert.h"
#include "core/atomic.h"
#include "core/octahedral.h"
#include "core/vec3f.h"

namespace ybi
{

#define YBI_PATH_GUIDING_CAMERA_VERTEX -1

struct PathGuidingTrainingSample
{
    Vec3f position;
    uint32_t octahedralDirection;
    float incidentRadiance;
    float pdf;
};

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

    YBI_CPU_GPU void InitializeIndirectSample(const Vec3f &samplePosition,
                                              const Vec3f &normal,
                                              const Vec3f &direction,
                                              int parentRecordIndex,
                                              const Vec3f &sampleThroughput,
                                              float samplePdf)
    {
        position = samplePosition;
        octahedralNormal = EncodeOctahedral(normal);
        octahedralDirection = EncodeOctahedral(direction);
        parentIndex = parentRecordIndex;
        references = 1u;
        accumulatedEmission = Vec3f(0.f);
        throughput = sampleThroughput;
        pdf = samplePdf;
    }

    YBI_CPU_GPU void InitializeDirectSample(const Vec3f &samplePosition,
                                            const Vec3f &normal,
                                            const Vec3f &direction,
                                            int parentRecordIndex,
                                            const Vec3f &sampleThroughput,
                                            float samplePdf,
                                            const Vec3f &sampleAccumulatedEmission)
    {
        InitializeIndirectSample(samplePosition, normal, direction, parentRecordIndex,
                                 sampleThroughput, samplePdf);
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
        if (recordIndex == YBI_PATH_GUIDING_CAMERA_VERTEX)
            return;
        RadianceRecord &record = GetRecord(recordIndex);
        if (AtomicSub(&record.references, 1u, AtomicOrder::AcqRel) == 1u)
        {
            BackpropagateRecursive(recordIndex);
        }
    }

    YBI_CPU_GPU bool WritePathGuidingTrainingSample(const RadianceRecord &record,
                                                    float incidentRadiance) const
    {
        const uint32_t sampleIndex = AtomicAdd(trainingSampleCount, 1u);
        if (sampleIndex >= trainingSamplesMax)
            return false;

        PathGuidingTrainingSample &sample = trainingSamples[sampleIndex];
        sample.position = record.position;
        sample.octahedralDirection = record.octahedralDirection;
        sample.incidentRadiance = incidentRadiance;
        sample.pdf = record.pdf;
        return true;
    }

    YBI_CPU_GPU void BackpropagateRecursive(int recordIndex) const
    {
        RadianceRecord &record = GetRecord(recordIndex);
        const Vec3f incidentRadiance(
            SafeDivide(record.accumulatedEmission.x, record.throughput.x),
            SafeDivide(record.accumulatedEmission.y, record.throughput.y),
            SafeDivide(record.accumulatedEmission.z, record.throughput.z));
        WritePathGuidingTrainingSample(record, incidentRadiance.MaxComponent());
        const int parentIndex = record.parentIndex;

        if (parentIndex != YBI_PATH_GUIDING_CAMERA_VERTEX)
        {
            RadianceRecord &parentRecord = GetRecord(parentIndex);
            parentRecord.AtomicAddEmission(record.accumulatedEmission);
            MarkDone(parentIndex);
        }
        Free(static_cast<uint32_t>(recordIndex));
    }
};

} // namespace ybi
