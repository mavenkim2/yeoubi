#include "render/integrator_ao.h"
#include "render/integrator_common.h"
#include "render/integrator_hit.h"
#include "render/integrator_path.h"
#include "render/integrator_primary.h"
#include "render/integrator_ray_differential.h"
#include "render/integrator_texture.h"
#include "render/launch_params.h"
#include "render/shading_core.h"
#include <assert.h>
#include <optix.h>
#include <optix_device.h>

using ybi::LaunchParams;
using ybi::render::integrator::HitInfo;
using ybi::render::integrator::UInt2;
using ybi::render::integrator::Vec3;
using ybi::render::integrator::Vec4;

extern "C"
{
    __constant__ LaunchParams params;
}

__device__ unsigned int g_try_sample_logged = 0u;

static __forceinline__ __device__ float3 Normalize3(const float3 &v)
{
    const float invLen = rsqrtf(v.x * v.x + v.y * v.y + v.z * v.z + 1e-20f);
    return float3{v.x * invLen, v.y * invLen, v.z * invLen};
}

static __forceinline__ __device__ float3 Cross3(const float3 &a, const float3 &b)
{
    return float3{a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}

static __forceinline__ __device__ Vec3 ToVec3(const float3 &v)
{
    return Vec3(v.x, v.y, v.z);
}

static __forceinline__ __device__ float3 ToFloat3(const Vec3 &v)
{
    return float3{v.x, v.y, v.z};
}

static __forceinline__ __device__ void
PackPointer(const void *ptr, unsigned int &payload0, unsigned int &payload1)
{
    const unsigned long long value = reinterpret_cast<unsigned long long>(ptr);
    payload0 = static_cast<unsigned int>(value & 0xffffffffu);
    payload1 = static_cast<unsigned int>((value >> 32u) & 0xffffffffu);
}

template <typename T>
static __forceinline__ __device__ T *UnpackPointer(unsigned int payload0, unsigned int payload1)
{
    const unsigned long long value = static_cast<unsigned long long>(payload0) |
                                     (static_cast<unsigned long long>(payload1) << 32u);
    return reinterpret_cast<T *>(value);
}

struct OptixState
{
    __device__ const LaunchParams &Params() const
    {
        return params;
    }

    __device__ UInt2 LaunchIndex() const
    {
        const uint3 idx = optixGetLaunchIndex();
        return {idx.x, idx.y};
    }

    __device__ bool SampleTexture2D(const LaunchParams::MaterialTextureRef &textureRef,
                                    float uu,
                                    float vv,
                                    Vec4 &outSample) const
    {
        const cudaTextureObject_t textureObject =
            static_cast<cudaTextureObject_t>(textureRef.textureObject);
        const float4 sample = tex2D<float4>(textureObject, uu, vv);
        outSample = {sample.x, sample.y, sample.z, sample.w};
        return true;
    }

    __device__ void RecordFeedbackKey(unsigned long long key) const
    {
        unsigned int *stats = reinterpret_cast<unsigned int *>(params.feedbackStats);
        unsigned long long *keys = reinterpret_cast<unsigned long long *>(params.feedbackKeys);
        const unsigned int index = atomicAdd(&stats[0], 1u);
        if (index < static_cast<unsigned int>(params.feedbackCapacity))
        {
            keys[index] = key;
        }
        else
        {
            atomicAdd(&stats[1], 1u);
        }
    }

    __device__ bool
    TraceOcclusion(const Vec3 &origin, const Vec3 &direction, float tMin, float tMax) const
    {
        unsigned int hit = 0x80000000u;
        const OptixTraversableHandle handle =
            static_cast<OptixTraversableHandle>(params.traversable);
        optixTrace(handle,
                   ToFloat3(origin),
                   ToFloat3(direction),
                   tMin,
                   tMax,
                   0.0f,
                   OptixVisibilityMask(255),
                   OPTIX_RAY_FLAG_DISABLE_CLOSESTHIT | OPTIX_RAY_FLAG_ENFORCE_ANYHIT |
                       OPTIX_RAY_FLAG_TERMINATE_ON_FIRST_HIT,
                   0,
                   1,
                   0,
                   hit);
        return hit != 0u;
    }

    __device__ bool TraceLightOcclusion(const Vec3 &origin,
                                        const Vec3 &direction,
                                        float tMin,
                                        float tMax,
                                        int /*lightIndex*/) const
    {
        return TraceOcclusion(origin, direction, tMin, tMax);
    }

    __device__ bool TraceClosest(
        const Vec3 &origin, const Vec3 &direction, float tMin, float tMax, HitInfo *outHit) const
    {
        if (!outHit)
        {
            return false;
        }
        *outHit = {};
        outHit->instanceId = -1;
        outHit->primitiveIndex = -1;

        unsigned int payload0 = 0u;
        unsigned int payload1 = 0u;
        PackPointer(outHit, payload0, payload1);

        const OptixTraversableHandle handle =
            static_cast<OptixTraversableHandle>(params.traversable);
        optixTrace(handle,
                   ToFloat3(origin),
                   ToFloat3(direction),
                   tMin,
                   tMax,
                   0.0f,
                   OptixVisibilityMask(255),
                   OPTIX_RAY_FLAG_NONE,
                   0,
                   1,
                   0,
                   payload0,
                   payload1);
        return outHit->instanceId >= 0;
    }

    __device__ void MaybeLogSampleSuccess() const
    {
        if (atomicCAS(&g_try_sample_logged, 0u, 1u) == 0u)
        {
            printf("TrySampleMaterialTexture succeeded\n");
        }
    }
};

static __forceinline__ __device__ unsigned int TraceColor(const float3 &origin,
                                                          const float3 &direction)
{
    unsigned int packedColor = 0;
    const OptixTraversableHandle handle = static_cast<OptixTraversableHandle>(params.traversable);
    optixTrace(handle,
               origin,
               direction,
               0.001f,
               1e20f,
               0.0f,
               OptixVisibilityMask(255),
               OPTIX_RAY_FLAG_NONE,
               0,
               1,
               0,
               packedColor);
    return packedColor;
}

static __forceinline__ __device__ bool
TryComputeTriangleNormal(int instanceId, int primitiveIndex, Vec3 &outNormal)
{
    if (params.instanceGeomRefs == 0ull || instanceId < 0 ||
        instanceId >= params.instanceGeomRefCount)
    {
        return false;
    }

    const LaunchParams::InstanceGeomRef *refs =
        reinterpret_cast<const LaunchParams::InstanceGeomRef *>(params.instanceGeomRefs);
    const LaunchParams::InstanceGeomRef ref = refs[instanceId];
    const int indexBase = primitiveIndex * 3;
    const float3 *positions = reinterpret_cast<const float3 *>(ref.positions);
    const int *indices = reinterpret_cast<const int *>(ref.indices);
    if (!positions || !indices || indexBase + 2 >= ref.numIndices)
    {
        return false;
    }

    const int i0 = indices[indexBase + 0];
    const int i1 = indices[indexBase + 1];
    const int i2 = indices[indexBase + 2];
    if (i0 < 0 || i0 >= ref.numPositions || i1 < 0 || i1 >= ref.numPositions || i2 < 0 ||
        i2 >= ref.numPositions)
    {
        return false;
    }

    const float3 p0 = positions[i0];
    const float3 p1 = positions[i1];
    const float3 p2 = positions[i2];
    const float3 edge01 = float3{p1.x - p0.x, p1.y - p0.y, p1.z - p0.z};
    const float3 edge02 = float3{p2.x - p0.x, p2.y - p0.y, p2.z - p0.z};
    const float3 localNormal = Normalize3(Cross3(edge01, edge02));
    const float3 worldNormal = Normalize3(optixTransformNormalFromObjectToWorldSpace(localNormal));
    outNormal = ToVec3(worldNormal);
    return true;
}

static __forceinline__ __device__ bool
TryComputeTriangleWorldPositions(int instanceId, int primitiveIndex, HitInfo *outHit)
{
    if (!outHit)
    {
        return false;
    }

    Vec3 p0 = {};
    Vec3 p1 = {};
    Vec3 p2 = {};
    if (!ybi::TryGetTriangleLocalPositions(params, instanceId, primitiveIndex, p0, p1, p2))
    {
        return false;
    }

    outHit->worldTri0 = ToVec3(optixTransformPointFromObjectToWorldSpace(ToFloat3(p0)));
    outHit->worldTri1 = ToVec3(optixTransformPointFromObjectToWorldSpace(ToFloat3(p1)));
    outHit->worldTri2 = ToVec3(optixTransformPointFromObjectToWorldSpace(ToFloat3(p2)));
    outHit->hasWorldTriangle = true;
    return true;
}

extern "C" __global__ void __raygen__primary()
{
    const uint3 launchIndex = optixGetLaunchIndex();
    const uint3 launchDims = optixGetLaunchDimensions();
    const unsigned int pixelIndex = launchIndex.y * launchDims.x + launchIndex.x;
    uchar4 *image = reinterpret_cast<uchar4 *>(params.image);
    if (!ybi::render::integrator::ShouldRenderLaunchPixel(params.singlePixelEnabled,
                                                          params.singlePixelX,
                                                          params.singlePixelY,
                                                          launchIndex.x,
                                                          launchIndex.y))
    {
        image[pixelIndex] = make_uchar4(0u, 0u, 0u, 255u);
        return;
    }
    const ybi::render::integrator::RayDifferential rayDiff =
        ybi::render::integrator::InitPerspectiveRayDifferential(
            params,
            static_cast<float>(launchIndex.x),
            static_cast<float>(launchIndex.y),
            launchDims.x,
            launchDims.y);
    const float3 centerDirection = ToFloat3(rayDiff.dir);
    const float3 origin = ToFloat3(rayDiff.origin);
    if (params.integrator == 3)
    {
        OptixState state = {};
        const unsigned int packedColor = ybi::render::integrator::IntegratorPathTrace(
            state, ToVec3(origin), ToVec3(centerDirection));
        image[pixelIndex] = make_uchar4((unsigned char)(packedColor & 255u),
                                        (unsigned char)((packedColor >> 8) & 255u),
                                        (unsigned char)((packedColor >> 16) & 255u),
                                        255u);
        return;
    }
    const unsigned int packedColor = TraceColor(origin, centerDirection);
    image[pixelIndex] = make_uchar4((unsigned char)(packedColor & 255u),
                                    (unsigned char)((packedColor >> 8) & 255u),
                                    (unsigned char)((packedColor >> 16) & 255u),
                                    255u);
}

extern "C" __global__ void __raygen__feedback()
{
    const uint3 launchIndex = optixGetLaunchIndex();
    const uint3 launchDims = optixGetLaunchDimensions();
    if (!ybi::render::integrator::ShouldRenderLaunchPixel(params.singlePixelEnabled,
                                                          params.singlePixelX,
                                                          params.singlePixelY,
                                                          launchIndex.x,
                                                          launchIndex.y))
    {
        return;
    }
    const ybi::render::integrator::RayDifferential rayDiff =
        ybi::render::integrator::InitPerspectiveRayDifferential(
            params,
            static_cast<float>(launchIndex.x),
            static_cast<float>(launchIndex.y),
            launchDims.x,
            launchDims.y);
    const float3 centerDirection = ToFloat3(rayDiff.dir);
    const float3 origin = ToFloat3(rayDiff.origin);
    TraceColor(origin, centerDirection);
}

extern "C" __global__ void __miss__primary()
{
    const unsigned int payload = optixGetPayload_0();
    if (payload == 0x80000000u)
    {
        optixSetPayload_0(0u);
        return;
    }

    if (params.integrator == 3)
    {
        return;
    }

    const float3 direction = optixGetWorldRayDirection();
    const Vec3 color = ybi::render::integrator::SkyColor(ToVec3(direction));
    optixSetPayload_0(ybi::render::PackRGB8(color.x, color.y, color.z));
}

extern "C" __global__ void __anyhit__primary()
{
    const unsigned int payload = optixGetPayload_0();
    if (payload == 0x80000000u)
    {
        optixSetPayload_0(1u);
        optixTerminateRay();
    }
}

extern "C" __global__ void __closesthit__primary()
{
    OptixState state;
    HitInfo hit = {};

    const float3 rayOrigin = optixGetWorldRayOrigin();
    const float3 rayDirection = Normalize3(optixGetWorldRayDirection());
    hit.rayOrigin = ToVec3(rayOrigin);
    hit.rayDir = ToVec3(rayDirection);
    hit.t = optixGetRayTmax();

    const unsigned int hitKind = optixGetHitKind();
    const bool isTriangle = (hitKind == OPTIX_HIT_KIND_TRIANGLE_FRONT_FACE ||
                             hitKind == OPTIX_HIT_KIND_TRIANGLE_BACK_FACE);

    if (isTriangle)
    {
        const float2 bary = optixGetTriangleBarycentrics();
        hit.barycentrics = Vec3(1.0f - bary.x - bary.y, bary.x, bary.y);
        hit.hasBarycentrics = true;
        hit.instanceId = static_cast<int>(optixGetInstanceId());
        hit.primitiveIndex = static_cast<int>(optixGetPrimitiveIndex());
        (void)TryComputeTriangleWorldPositions(hit.instanceId, hit.primitiveIndex, &hit);
        (void)ybi::TryComputeTriangleSurfacePartials(params, &hit);

        Vec3 geomNormal(-rayDirection.x, -rayDirection.y, -rayDirection.z);
        if (TryComputeTriangleNormal(hit.instanceId, hit.primitiveIndex, geomNormal))
        {
            hit.geomNormal = geomNormal;
            hit.hasGeomNormal = true;
        }

        hit.hasShadingNormal = ybi::ComputeTriangleShadingNormal(params,
                                                                 hit.barycentrics.y,
                                                                 hit.barycentrics.z,
                                                                 hit.primitiveIndex,
                                                                 hit.instanceId,
                                                                 hit.shadingNormal);
        if (hit.hasShadingNormal)
        {
            hit.shadingNormal = Normalize(ToVec3(
                optixTransformNormalFromObjectToWorldSpace(
                    {hit.shadingNormal.x, hit.shadingNormal.y, hit.shadingNormal.z})));
        }

        if (hit.hasGeomNormal)
        {
            const ybi::render::integrator::RayDifferential *rayDiff = nullptr;
            ybi::render::integrator::RayDifferential primaryRayDiff = {};
            if (params.integrator != 3)
            {
                const uint3 launchIndex = optixGetLaunchIndex();
                const uint3 launchDims = optixGetLaunchDimensions();
                primaryRayDiff = ybi::render::integrator::InitPerspectiveRayDifferential(
                    params,
                    static_cast<float>(launchIndex.x),
                    static_cast<float>(launchIndex.y),
                    launchDims.x,
                    launchDims.y);
                rayDiff = &primaryRayDiff;
            }
            (void)ybi::TryComputeTriangleHitDifferentials(params, rayDiff, &hit);
        }
    }

    if (params.integrator == 3)
    {
        HitInfo *outHit = UnpackPointer<HitInfo>(optixGetPayload_0(), optixGetPayload_1());
        if (outHit)
        {
            *outHit = hit;
        }
        return;
    }

    if (params.integrator == 2)
    {
        if (hit.hasBarycentrics)
        {
            ybi::render::integrator::IntegratorFeedbackOnly(state, hit);
        }
        optixSetPayload_0(0u);
        return;
    }

    if (params.integrator == 1)
    {
        Vec3 geomNormal(-rayDirection.x, -rayDirection.y, -rayDirection.z);
        bool hasGeom = false;
        if (isTriangle && hit.hasBarycentrics && !hit.hasGeomNormal)
        {
            hasGeom = TryComputeTriangleNormal(hit.instanceId, hit.primitiveIndex, geomNormal);
            hit.geomNormal = geomNormal;
            hit.hasGeomNormal = hasGeom;
        }
        else
        {
            geomNormal = hit.geomNormal;
            hasGeom = hit.hasGeomNormal;
        }
        hit.geomNormal = geomNormal;
        hit.hasGeomNormal = hasGeom;
        const unsigned int packed = ybi::render::integrator::IntegratorAO(state, hit);
        optixSetPayload_0(packed);
        return;
    }

    if (!hit.hasBarycentrics)
    {
        optixSetPayload_0(ybi::render::PackRGB8(0.7f, 0.7f, 0.7f));
        return;
    }

    const unsigned int packed = ybi::render::integrator::IntegratorPrimaryDiffuse(state, hit);
    optixSetPayload_0(packed);
}
