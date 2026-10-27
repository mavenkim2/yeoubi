#pragma once

#include "core/math.h"
#include "device/cuda_compat.h"
#include "kernel/backend/dispatch.h"
#if defined(YBI_ENABLE_PATH_GUIDING)
#include "kernel/guiding/api.h"
#endif
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
    const Vec3f position = sampleState->ReadPosition(sampleIndex);
    RNG rng = sampleState->ReadRNG(sampleIndex);
    BsdfSample sample;

#ifdef YBI_ENABLE_PATH_GUIDING
    auto sampleBsdf = [&](BsdfSample &sample) {
        const bool success = DeviceDispatch<bool>(material.sampleCallableIndex,
                                                  &params,
                                                  &closureInfo,
                                                  viewDirection,
                                                  &sample,
                                                  Vec3f(rng.Uniform(), rng.Uniform2D()));
        return success;
    };
    auto evalBsdf = [&](const Vec3f &direction) {
        BsdfSample sample{};
        BsdfEvaluationProperties properties;

        const Vec3f bsdfCos = DeviceDispatch<Vec3f>(material.evaluateCallableIndex,
                                                    &params,
                                                    &closureInfo,
                                                    viewDirection,
                                                    direction,
                                                    &properties);
        sample.direction = direction;
        sample.pdf = DeviceDispatch<float>(
            material.pdfCallableIndex, &params, &closureInfo, viewDirection, direction);
        if (sample.pdf > 0.f)
            sample.weight = bsdfCos / sample.pdf;
        // The current material evaluation callable supports continuous reflection.
        sample.eta = 1.f;
        sample.flags = properties.flags;
        sample.sampledRoughness = properties.maxRoughness;

        return sample;
    };

    const bool success =
        SampleGuiding(sample, params.pathGuiding, position, sampleBsdf, evalBsdf, rng.Uniform3D());
#else
    const bool success = DeviceDispatch<bool>(material.sampleCallableIndex,
                                              &params,
                                              &closureInfo,
                                              viewDirection,
                                              &sample,
                                              Vec3f(rng.Uniform(), rng.Uniform2D()));
#endif

    sampleState->WriteRNG(sampleIndex, rng);
    if (!success)
    {
        TerminatePath(params, sampleState, sampleIndex);
        return;
    }

    const Vec3f nextDirection = sample.direction;
    Vec3f throughput = sampleState->ReadThroughput(sampleIndex) * sample.weight;
    const float etaScale = sampleState->ReadEtaScale(sampleIndex) * Sqr(sample.eta);
#if defined(YBI_ENABLE_PATH_GUIDING)
    const double russianRouletteCorrectionFactor = UpdateRussianRouletteCorrectionFactor(
        sampleState->ReadRussianRouletteCorrectionFactor(sampleIndex), sample);
    sampleState->WriteRussianRouletteCorrectionFactor(sampleIndex, russianRouletteCorrectionFactor);
    const float survivalProbability = float(Min(
        double(throughput.MaxComponent()) * etaScale * russianRouletteCorrectionFactor, 1.0));
#else
    const float survivalProbability = (throughput * etaScale).MaxComponent();
#endif
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

    Vec3f geometricNormal = sampleState->ReadGeometricNormal(sampleIndex);
    if (Dot(nextDirection, geometricNormal) < 0.f)
        geometricNormal = -geometricNormal;

#ifdef YBI_ENABLE_PATH_GUIDING
    const int recordIndex = sampleState->ReadRadianceRecordIndex(sampleIndex);
    const int childRecordIndex = params.radianceRecorder.RecordIndirectSample(
        position, sampleState->ReadShadingNormal(sampleIndex), sample, recordIndex, throughput);
    sampleState->WriteRadianceRecordIndex(sampleIndex, childRecordIndex);
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
