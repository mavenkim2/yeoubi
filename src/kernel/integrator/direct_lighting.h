#pragma once

#include "core/math.h"
#include "device/cuda_compat.h"
#include "kernel/backend/dispatch.h"
#if defined(YBI_ENABLE_PATH_GUIDING)
#include "kernel/guiding/api.h"
#endif
#include "kernel/integrator/mis.h"
#include "kernel/launch/kernel_params.h"
#include "kernel/light/sample.h"
#include "kernel/light/visibility.h"
#include "kernel/state/types.h"

namespace ybi
{

template <typename SampleState>
YBI_GPU void
DirectLighting(const KernelParams &params, const SampleState *sampleState, uint32_t sampleIndex)
{
    sampleState->WriteKernelType(sampleIndex, KERNEL_INDIRECT_LIGHTING);

    const DeviceScene &scene = params.scene;
    if (GetLightCount(scene) == 0)
        return;

    const uint32_t materialIndex = sampleState->ReadMaterialIndex(sampleIndex);

    const MaterialHandle &handle = scene.materialHandles[materialIndex];
    if (handle.type != MaterialHandleType::MaterialX)
        return;

    const DeviceMaterial3 &material = scene.deviceMaterial3s[handle.index];
    const Vec3f position = sampleState->ReadPosition(sampleIndex);
    const Vec3f geometricNormal = sampleState->ReadGeometricNormal(sampleIndex);
    const Vec3f viewDirection = -sampleState->ReadRayDirection(sampleIndex);
    const Vec3f throughput = sampleState->ReadThroughput(sampleIndex);
    Vec3f radiance = sampleState->ReadRadiance(sampleIndex);
    RNG rng = sampleState->ReadRNG(sampleIndex);
    const ClosureEvaluationInfo closureInfo = sampleState->ReadClosureEvaluationInfo(sampleIndex);

    LightSample lightSample;
    if (SampleLight(scene, position, Vec3f(rng.Uniform(), rng.Uniform2D()), lightSample))
    {
        BsdfEvaluationProperties *evaluationProperties = nullptr;
#if defined(YBI_ENABLE_PATH_GUIDING)
        BsdfEvaluationProperties properties;
        evaluationProperties = &properties;
#endif
        const Vec3f sampledLightBsdfWeight = DeviceDispatch<Vec3f>(
            material.evaluateCallableIndex,
            &params,
            &closureInfo,
            viewDirection,
            lightSample.direction,
            evaluationProperties);
        const float sampledLightBsdfPdf = DeviceDispatch<float>(material.pdfCallableIndex,
                                                                &params,
                                                                &closureInfo,
                                                                viewDirection,
                                                                lightSample.direction);
        if (sampledLightBsdfWeight.MaxComponent() != 0.f)
        {
            Vec3f contribution(0.f);
            if (IsLightSampleVisible(params.bvhHandle, lightSample, position, geometricNormal))
            {
                float continuationPdf = sampledLightBsdfPdf;
#if defined(YBI_ENABLE_PATH_GUIDING)
                continuationPdf = GuidingProposalPDF(
                    params.pathGuiding, position, lightSample.direction, sampledLightBsdfPdf);
#endif
                const float misWeight = PowerHeuristic(lightSample.pdf, continuationPdf);
                contribution = throughput * lightSample.radiance * sampledLightBsdfWeight *
                               misWeight / lightSample.pdf;
                radiance += contribution;
            }
#if defined(YBI_ENABLE_PATH_GUIDING)
            const int parentRecordIndex = sampleState->ReadRadianceRecordIndex(sampleIndex);
            params.radianceRecorder.RecordDirectSample(position,
                                                       sampleState->ReadShadingNormal(sampleIndex),
                                                       lightSample.direction,
                                                       parentRecordIndex,
                                                       throughput * sampledLightBsdfWeight /
                                                           lightSample.pdf,
                                                       lightSample.pdf,
                                                       lightSample.distance,
                                                       contribution,
                                                       properties);
#endif
        }
    }

    sampleState->WriteRNG(sampleIndex, rng);
    sampleState->WriteRadiance(sampleIndex, radiance);
}

} // namespace ybi
