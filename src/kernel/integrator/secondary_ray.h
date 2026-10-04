#pragma once

#include "kernel/geometry/triangle_info.h"
#include "kernel/launch/kernel_params.h"
#include "kernel/raytracing/ray_cone.h"
#include "kernel/raytracing/trace.h"
#include "kernel/state/types.h"

namespace ybi
{

template <typename SampleState>
YBI_GPU void
SecondaryRay(const KernelParams &params, const SampleState *sampleState, uint32_t sampleIndex)
{
    Vec3f rayOrigin = sampleState->ReadRayOrigin(sampleIndex);
    Vec3f rayDirection = sampleState->ReadRayDirection(sampleIndex);
    RayCone rayCone = sampleState->ReadRayCone(sampleIndex);

    const IntersectionInfo info =
        TraceRay(params.bvhHandle, rayOrigin, rayDirection, 0.f, 1.0e30f);

    if (info.type != INTERSECTION_TYPE_NONE)
    {
        const Object &object = params.scene.objects[info.id];
        TriangleInfo triangleInfo;
        GetTriangleInfo(params.scene, object, info, triangleInfo);
        triangleInfo.TransformToWorld(info.objectToWorld);

        TriangleSurfaceInfo surfaceInfo;
        CalculateTriangleSurfaceInfo(triangleInfo, info.bary, surfaceInfo);
        UpdateSecondaryRayCone(rayCone, rayDirection, info.t, triangleInfo, surfaceInfo);

        sampleState->WriteRayCone(sampleIndex, rayCone);
        sampleState->WriteIntersection(sampleIndex, info);
        sampleState->WriteTriangleInfo(sampleIndex, surfaceInfo);
        sampleState->WriteKernelType(sampleIndex, KERNEL_SHADE_HIT);
    }
    else
    {
        Vec3f radiance = sampleState->ReadRadiance(sampleIndex);
        Vec3f throughput = sampleState->ReadThroughput(sampleIndex);

        // TODO: stop hardcoding this :)
        Vec3f environment(0.5f);
        radiance += throughput * environment;

        Vec3f *image = params.kernelOutput;
        const uint32_t pixelIndex = sampleState->ReadPixelIndex(sampleIndex);
        image[pixelIndex] = radiance;
        sampleState->WriteRadiance(sampleIndex, radiance);
        sampleState->WriteKernelType(sampleIndex, KERNEL_FINISHED);
    }
}

} // namespace ybi
