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
    Vec3f nextDirection;

    const bool success = DeviceDispatch<bool>(material.sampleCallableIndex,
                                              &params,
                                              &closureInfo,
                                              viewDirection,
                                              &nextDirection,
                                              Vec3f(rng.Uniform(), rng.Uniform2D()));
    sampleState->WriteRNG(sampleIndex, rng);
    if (!success)
    {
        TerminatePath(params, sampleState, sampleIndex);
        return;
    }

    const float pdf = DeviceDispatch<float>(
        material.pdfCallableIndex, &params, &closureInfo, viewDirection, nextDirection);
    if (!(pdf > 0.f))
    {
        TerminatePath(params, sampleState, sampleIndex);
        return;
    }

    const Vec3f weight =
        DeviceDispatch<Vec3f>(
            material.evaluateCallableIndex, &params, &closureInfo, viewDirection, nextDirection) /
        pdf;
    Vec3f throughput = sampleState->ReadThroughput(sampleIndex) * weight;

    // Match the roulette policy used by the unidirectional integrator.
    const float survivalProbability = throughput.MaxComponent();
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
                                   pdf);

    sampleState->WriteRadianceRecordIndex(sampleIndex, childRecordIndex);
    params.radianceRecorder.IncrementCount(recordIndex);
    params.radianceRecorder.MarkDone(recordIndex);
#endif

    sampleState->WriteThroughput(sampleIndex, throughput);
    sampleState->WritePreviousSurfacePosition(sampleIndex, position);
    sampleState->WritePreviousBsdfPdf(sampleIndex, pdf);
    sampleState->WriteRayOrigin(sampleIndex, OffsetRayOrigin(position, geometricNormal));
    sampleState->WriteRayDirection(sampleIndex, nextDirection);
    sampleState->WriteDepth(sampleIndex, depth + 1);
    sampleState->WriteKernelType(sampleIndex, KERNEL_INTERSECT);
}

} // namespace ybi
