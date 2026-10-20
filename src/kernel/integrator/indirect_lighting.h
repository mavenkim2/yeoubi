#pragma once

#include "core/math.h"
#include "device/cuda_compat.h"
#include "kernel/backend/dispatch.h"
#include "kernel/integrator/terminate_path.h"
#include "kernel/launch/kernel_params.h"
#include "kernel/raytracing/ray_offset.h"
#include "kernel/state/types.h"
#include <cassert>

namespace ybi
{

template <typename SampleState>
YBI_GPU void
IndirectLighting(const KernelParams &params, const SampleState *sampleState, uint32_t sampleIndex)
{
    const uint32_t depth = sampleState->ReadDepth(sampleIndex);

    const DeviceScene &scene = params.scene;
    const uint32_t materialIndex = sampleState->ReadMaterialIndex(sampleIndex);

    const MaterialHandle &handle = scene.materialHandles[materialIndex];

    assert(handle.type == MaterialHandleType::MaterialX);

    const DeviceMaterial3 &material = scene.deviceMaterial3s[handle.index];
    const ClosureEvaluationInfo closureInfo = sampleState->ReadClosureEvaluationInfo(sampleIndex);
    const Vec3f viewDirection = -sampleState->ReadRayDirection(sampleIndex);
    RNG rng = sampleState->ReadRNG(sampleIndex);
    BsdfSample sample;

    const bool success = DeviceDispatch<bool>(material.sampleCallableIndex,
                                              &params,
                                              &closureInfo,
                                              viewDirection,
                                              &sample,
                                              Vec3f(rng.Uniform(), rng.Uniform2D()));
    sampleState->WriteRNG(sampleIndex, rng);
    if (!success)
    {
        TerminatePath(params, sampleState, sampleIndex);
        return;
    }

    const Vec3f nextDirection = sample.direction;
    Vec3f throughput = sampleState->ReadThroughput(sampleIndex) * sample.weight;
    const float etaScale = sampleState->ReadEtaScale(sampleIndex) * Sqr(sample.eta);
    const float survivalProbability = (throughput * etaScale).MaxComponent();
    if (survivalProbability < 1.f && depth > 1)
    {
        const float rouletteSample = rng.Uniform();
        sampleState->WriteRNG(sampleIndex, rng);
        if (rouletteSample >= survivalProbability)
        {
            sampleState->WriteThroughput(sampleIndex, throughput);
            TerminatePath(params, sampleState, sampleIndex);
            return;
        }
        throughput /= survivalProbability;
    }

    const Vec3f position = sampleState->ReadPosition(sampleIndex);
    Vec3f geometricNormal = sampleState->ReadGeometricNormal(sampleIndex);
    if (Dot(nextDirection, geometricNormal) < 0.f)
        geometricNormal = -geometricNormal;

#ifdef YBI_USE_PATH_GUIDING
    int recordIndex = sampleState->ReadRadianceRecordIndex(sampleIndex);
    int childRecordIndex;
    const bool recordAllocated = params.radianceRecorder.AllocateRecord(childRecordIndex);
    assert(recordAllocated);

    RadianceRecord &child = params.radianceRecorder.records[childRecordIndex];
    child.InitializeIndirectSample(position,
                                   sampleState->ReadShadingNormal(sampleIndex),
                                   nextDirection,
                                   recordIndex,
                                   throughput,
                                   sample.pdf);

    sampleState->WriteRadianceRecordIndex(sampleIndex, childRecordIndex);
    params.radianceRecorder.IncrementCount(recordIndex);
    params.radianceRecorder.MarkDone(recordIndex);
#endif

    sampleState->WriteThroughput(sampleIndex, throughput);
    sampleState->WriteEtaScale(sampleIndex, etaScale);
    sampleState->WritePreviousSurfacePosition(sampleIndex, position);
    sampleState->WritePreviousBsdfPdf(sampleIndex, sample.pdf);
    sampleState->WritePreviousBsdfFlags(sampleIndex, sample.flags);
    sampleState->WriteRayOrigin(sampleIndex, OffsetRayOrigin(position, geometricNormal));
    sampleState->WriteRayDirection(sampleIndex, nextDirection);
    sampleState->WriteDepth(sampleIndex, depth + 1);
    sampleState->WriteKernelType(sampleIndex, KERNEL_INTERSECT);
}

} // namespace ybi
