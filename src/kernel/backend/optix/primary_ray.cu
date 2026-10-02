#include "kernel/geometry/triangle_info.h"
#include "kernel/integrator/ray_cone.h"
#include "kernel/launch/wavefront_params.h"
#include "kernel/sampling/rng.h"
#include "kernel/state/types.h"
#include <optix_device.h>

namespace ybi
{

extern "C"
{
    __constant__ WavefrontParams params;
}

YBI_GPU void PrimaryRayKernel()
{
    const uint3 launchIndex = optixGetLaunchIndex();
    const uint3 launchDimensions = optixGetLaunchDimensions();

    const uint32_t rasterX = launchIndex.x + params.pixelOffsetX;
    const uint32_t rasterY = launchIndex.y + params.pixelOffsetY;
    const OptixTraversableHandle bvhHandle = static_cast<OptixTraversableHandle>(params.bvhHandle);

    RNG rng = RNG::Init(Vec2u(rasterX, rasterY), params.sampleOffset);
    const float pixelX = float(rasterX) + rng.Uniform();
    const float pixelY = float(rasterY) + rng.Uniform();
    Vec4f cameraPos = params.cameraFromRaster * Vec4f(Vec3f(pixelX, pixelY, 0.f), 1.f);
    cameraPos.xyz /= cameraPos.w;

    Vec3f rayOrigin = params.renderFromCamera.TransformPoint(Vec3f(0.f));
    Vec3f rayDirection = Normalize(params.renderFromCamera.TransformVector(cameraPos.xyz));

    RayCone rayCone{};
    rayCone.spreadAngle = atanf(2.f * tanf(0.5f * params.fovRadians) / float(params.imageHeight));

    IntersectionInfo info{};
    const uint64_t ptr = reinterpret_cast<uint64_t>(&info);
    uint32_t low = static_cast<uint32_t>(ptr);
    uint32_t high = static_cast<uint32_t>(ptr >> 32u);

    optixTrace(bvhHandle,
               make_float3(rayOrigin.x, rayOrigin.y, rayOrigin.z),
               make_float3(rayDirection.x, rayDirection.y, rayDirection.z),
               0.f,
               1.0e30f,
               0.f,
               0xff,
               OPTIX_RAY_FLAG_NONE,
               0,
               0,
               0,
               low,
               high);

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
        params.sampleState.WritePixelX(sampleIndex, rasterX);
        params.sampleState.WritePixelY(sampleIndex, rasterY);
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
        const uint32_t pixelIndex = launchIndex.x + launchIndex.y * launchDimensions.x;
        image[pixelIndex] = radiance;
    }
}

extern "C" __global__ void __miss__kernel()
{
    uint32_t low = optixGetPayload_0();
    uint32_t high = optixGetPayload_1();
    uint64_t ptr = low | (uint64_t(high) << 32u);

    IntersectionInfo *isect = reinterpret_cast<IntersectionInfo *>(ptr);
    isect->type = INTERSECTION_TYPE_NONE;
}

extern "C" __global__ void __closesthit__kernel()
{
    uint32_t low = optixGetPayload_0();
    uint32_t high = optixGetPayload_1();
    uint64_t ptr = low | (uint64_t(high) << 32u);

    IntersectionInfo *isect = reinterpret_cast<IntersectionInfo *>(ptr);

    if (optixIsTriangleHit())
    {
        const float2 bary = optixGetTriangleBarycentrics();
        float objectToWorld[12];
        // float worldToObject[12];
        optixGetObjectToWorldTransformMatrix(objectToWorld);
        // optixGetWorldToObjectTransformMatrix(worldToObject);

        isect->t = optixGetRayTmax();
        isect->id = optixGetInstanceId();
        isect->primIndex = optixGetPrimitiveIndex();
        isect->clusterId = ~0u;
#if OPTIX_VERSION >= 90000
        isect->clusterId = optixGetClusterId();
#endif
        isect->bary = Vec2f(bary.x, bary.y);
        isect->objectToWorld = AffineTransform(objectToWorld);
        // isect->worldToObject = AffineTransform(worldToObject);
        isect->type = INTERSECTION_TYPE_TRIANGLE;
    }
    else
    {
    }
}

extern "C" __global__ void __raygen__wavefront()
{
    PrimaryRayKernel();
}

} // namespace ybi
