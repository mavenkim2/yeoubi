#pragma once

#include "kernel/geometry/triangle_info.h"
#include "kernel/integrator/mis.h"
#include "kernel/integrator/terminate_path.h"
#include "kernel/launch/kernel_params.h"
#include "kernel/light/envmap.h"
#include "kernel/light/sample.h"
#include "kernel/raytracing/ray_cone.h"
#include "kernel/raytracing/ray_offset.h"
#include "kernel/raytracing/trace.h"
#include "kernel/state/types.h"

namespace ybi
{

template <typename SampleState>
YBI_GPU void
Intersect(const KernelParams &params, const SampleState *sampleState, uint32_t sampleIndex)
{
    Vec3f rayOrigin = sampleState->ReadRayOrigin(sampleIndex);
    const Vec3f rayDirection = sampleState->ReadRayDirection(sampleIndex);
    RayCone rayCone = sampleState->ReadRayCone(sampleIndex);
    Vec3f radiance = sampleState->ReadRadiance(sampleIndex);
    const Vec3f throughput = sampleState->ReadThroughput(sampleIndex);
    const uint32_t depth = sampleState->ReadDepth(sampleIndex);

    const IntersectionInfo info =
        TraceRay(params.bvhHandle, rayOrigin, rayDirection, 0.f, 1.0e30f);

    if (info.type == INTERSECTION_TYPE_NONE)
    {
        Vec3f environment(0.f);
        for (uint32_t lightIndex = 0; lightIndex < params.scene.numEnvmapLights; lightIndex++)
            environment +=
                EvaluateEnvmapLight(params.scene.envmapLights[lightIndex], rayDirection);
        const Vec3f contribution = throughput * environment;
        radiance += contribution;
#if defined(YBI_USE_PATH_GUIDING)
        const int recordIndex = sampleState->ReadRadianceRecordIndex(sampleIndex);
        if (recordIndex != YBI_PATH_GUIDING_CAMERA_VERTEX)
            params.radianceRecorder.GetRecord(recordIndex).AtomicAddEmission(contribution);
#endif
        TerminatePath(params, sampleState, sampleIndex, radiance);
        return;
    }
    else
    {
        const Object &object = params.scene.objects[info.id];
        if (object.lightIndex >= 0)
        {
            const DeviceLight &light = params.scene.deviceLights[object.lightIndex];
            const Vec3f lightHitPosition = rayOrigin + info.t * rayDirection;

            // Rect lights emit from local -Z; back-face hits pass through.
            const Vec3f localRayDirection = light.localFromWorld.TransformVector(rayDirection);
            if (localRayDirection.z > 0.f)
            {
                float misWeight = 1.f;
                const float previousBsdfPdf = sampleState->ReadPreviousBsdfPdf(sampleIndex);
                if (previousBsdfPdf > 0.f)
                {
                    const float lightPdf =
                        PDFLight(params.scene,
                                 static_cast<uint32_t>(object.lightIndex),
                                 sampleState->ReadPreviousSurfacePosition(sampleIndex),
                                 lightHitPosition);
                    misWeight = PowerHeuristic(previousBsdfPdf, lightPdf);
                }
                radiance += throughput * light.radiance * misWeight;
            }

            const Vec3f lightTangentX = info.objectToWorld.TransformVector(Vec3f(1.f, 0.f, 0.f));
            const Vec3f lightTangentY = info.objectToWorld.TransformVector(Vec3f(0.f, 1.f, 0.f));
            Vec3f lightNormal = Normalize(Cross(lightTangentX, lightTangentY));
            if (Dot(lightNormal, rayDirection) < 0.f)
            {
                lightNormal = -lightNormal;
            }
            rayOrigin = OffsetRayOrigin(lightHitPosition, lightNormal);
            sampleState->WriteKernelType(sampleIndex, KERNEL_INTERSECT);
        }
        else if (depth >= static_cast<uint32_t>(Max(0, params.maxDepth)))
        {
            TerminatePath(params, sampleState, sampleIndex);
            return;
        }
        else
        {
            TriangleInfo triangleInfo;
            GetTriangleInfo(params.scene, object, info, triangleInfo);
            CalculateTriangleInfo(triangleInfo, info.bary);
            triangleInfo.TransformToWorld(info.objectToWorld);

            TriangleSurfaceInfo surfaceInfo;
            CalculateTriangleSurfaceInfo(triangleInfo, surfaceInfo);
            if (depth == 0)
            {
                UpdatePrimaryRayCone(rayCone, rayDirection, info.t, triangleInfo, surfaceInfo);
            }
            else
            {
                UpdateSecondaryRayCone(rayCone, rayDirection, info.t, triangleInfo, surfaceInfo);
            }

            sampleState->WriteIntersection(sampleIndex, info);
            sampleState->WriteTriangleInfo(sampleIndex, surfaceInfo);
            sampleState->WriteKernelType(sampleIndex, KERNEL_SHADE_HIT);
        }
    }

    sampleState->WriteRayOrigin(sampleIndex, rayOrigin);
    sampleState->WriteRayCone(sampleIndex, rayCone);
    sampleState->WriteRadiance(sampleIndex, radiance);
}

} // namespace ybi
