#pragma once

#include "core/vec2u.h"
#include "kernel/geometry/triangle_info.h"
#include "kernel/launch/wavefront_params.h"
#include "kernel/raytracing/ray_cone.h"
#include "kernel/raytracing/trace.h"
#include "kernel/sampling/rng.h"
#include "kernel/state/types.h"

namespace ybi
{

YBI_GPU void
PrimaryRay(const WavefrontParams &params, const Vec2u &launchIndex, const Vec2u &launchDimensions)
{
    const uint32_t rasterX = launchIndex.x + params.pixelOffsetX;
    const uint32_t rasterY = launchIndex.y + params.pixelOffsetY;
    const uint32_t pixelIndex = launchIndex.x + launchIndex.y * launchDimensions.x;

    RNG rng = RNG::Init(Vec2u(rasterX, rasterY), params.sampleOffset);
    const float pixelX = float(rasterX) + rng.Uniform();
    const float pixelY = float(rasterY) + rng.Uniform();
    Vec4f cameraPos = params.cameraFromRaster * Vec4f(Vec3f(pixelX, pixelY, 0.f), 1.f);
    cameraPos.xyz /= cameraPos.w;

    Vec3f rayOrigin = params.renderFromCamera.TransformPoint(Vec3f(0.f));
    Vec3f rayDirection = Normalize(params.renderFromCamera.TransformVector(cameraPos.xyz));

    RayCone rayCone{};
    rayCone.spreadAngle = atanf(2.f * tanf(0.5f * params.fovRadians) / float(params.imageHeight));

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
        UpdatePrimaryRayCone(rayCone, rayDirection, info.t, triangleInfo, surfaceInfo);

        uint32_t sampleIndex = atomicAdd(params.numSampleStates, 1u);
        params.sampleState.WritePixelIndex(sampleIndex, pixelIndex);
        params.sampleState.WriteRNG(sampleIndex, rng);
        params.sampleState.WriteRayCone(sampleIndex, rayCone);
        params.sampleState.WriteRadiance(sampleIndex, Vec3f(0.f));
        params.sampleState.WriteThroughput(sampleIndex, Vec3f(1.f));
        params.sampleState.WriteRayOrigin(sampleIndex, rayOrigin);
        params.sampleState.WriteRayDirection(sampleIndex, rayDirection);
        params.sampleState.WriteDepth(sampleIndex, 0u);
        params.sampleState.WriteIntersection(sampleIndex, info);
        params.sampleState.WriteTriangleInfo(sampleIndex, surfaceInfo);
        params.sampleState.WriteKernelType(sampleIndex, KERNEL_SHADE_HIT);
    }
    else
    {
        // TODO: stop hardcoding this :)
        Vec3f radiance(0.5f);
        Vec3f *image = params.kernelOutput;
        image[pixelIndex] = radiance;
    }
}

} // namespace ybi
