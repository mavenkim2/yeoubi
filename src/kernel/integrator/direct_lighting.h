#pragma once

#include "core/math.h"
#include "device/cuda_compat.h"
#include "kernel/backend/dispatch.h"
#include "kernel/integrator/mis.h"
#include "kernel/launch/kernel_params.h"
#include "kernel/light/sample.h"
#include "kernel/raytracing/ray_offset.h"
#include "kernel/raytracing/trace.h"
#include "kernel/state/types.h"

namespace ybi
{

template <typename SampleState>
YBI_GPU void
DirectLighting(const KernelParams &params, const SampleState *sampleState, uint32_t sampleIndex)
{
    sampleState->WriteKernelType(sampleIndex, KERNEL_INDIRECT_LIGHTING);

    const DeviceScene &scene = params.scene;
    if (scene.numLights == 0)
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

    int sampledLightIndex;
    Vec3f sampledLightPosition;
    float lightPdf;

    const Vec3f sampledLightDirection = SampleLight(scene,
                                                    position,
                                                    Vec3f(rng.Uniform(), rng.Uniform2D()),
                                                    sampledLightIndex,
                                                    sampledLightPosition,
                                                    lightPdf);

    if (lightPdf != 0.f)
    {
        const DeviceLight &sampledLight = scene.deviceLights[sampledLightIndex];

        const Vec3f sampledLightBsdfWeight = DeviceDispatch<Vec3f>(material.evaluateCallableIndex,
                                                                   &params,
                                                                   &closureInfo,
                                                                   viewDirection,
                                                                   sampledLightDirection);
        const float sampledLightBsdfPdf = DeviceDispatch<float>(material.pdfCallableIndex,
                                                                &params,
                                                                &closureInfo,
                                                                viewDirection,
                                                                sampledLightDirection);
        if (sampledLightBsdfWeight.MaxComponent() != 0.f)
        {
            const Vec3f shadowNormal = Dot(sampledLightDirection, geometricNormal) >= 0.f
                                           ? geometricNormal
                                           : -geometricNormal;
            const Vec3f shadowOrigin = OffsetRayOrigin(position, shadowNormal);
            const Vec3f shadowToLight = sampledLightPosition - shadowOrigin;
            const float shadowDistance = Length(shadowToLight);
            if (shadowDistance > 0.f)
            {
                const Vec3f shadowDirection = shadowToLight / shadowDistance;
                Vec3f contribution(0.f);
                if (TraceVisibility(params.bvhHandle,
                                    shadowOrigin,
                                    shadowDirection,
                                    shadowDistance * (1.f - 1e-4f)))
                {
                    const float misWeight = PowerHeuristic(lightPdf, sampledLightBsdfPdf);
                    contribution = throughput * sampledLight.radiance * sampledLightBsdfWeight *
                                   misWeight / lightPdf;
                    radiance += contribution;
                }
#if defined(YBI_USE_PATH_GUIDING)
                // TODO: do i want to create a sample here? could just add to previous node
                const int parentRecordIndex = sampleState->ReadRadianceRecordIndex(sampleIndex);
                int childRecordIndex;
                if (params.radianceRecorder.AllocateRecord(childRecordIndex))
                {
                    RadianceRecord &child = params.radianceRecorder.GetRecord(childRecordIndex);
                    child.InitializeDirectSample(position,
                                                 sampleState->ReadShadingNormal(sampleIndex),
                                                 sampledLightDirection,
                                                 parentRecordIndex,
                                                 throughput * sampledLightBsdfWeight / lightPdf,
                                                 lightPdf,
                                                 contribution);
                    params.radianceRecorder.IncrementCount(parentRecordIndex);
                    params.radianceRecorder.MarkDone(childRecordIndex);
                }
                else if (parentRecordIndex != YBI_PATH_GUIDING_CAMERA_VERTEX)
                {
                    params.radianceRecorder.GetRecord(parentRecordIndex)
                        .AtomicAddEmission(contribution);
                }
#endif
            }
        }
    }

    sampleState->WriteRNG(sampleIndex, rng);
    sampleState->WriteRadiance(sampleIndex, radiance);
}

} // namespace ybi
