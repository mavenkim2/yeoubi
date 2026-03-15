#include "render/cpu_kernels.h"

#if defined(WITH_EMBREE)

#include "device/cpu_device.h"
#include "render/integrator_ao.h"
#include "render/integrator_hit.h"
#include "render/integrator_path.h"
#include "render/integrator_primary.h"
#include "render/integrator_ray_differential.h"
#include "render/launch_params.h"
#include "render/shading_core.h"
#include "texture/virtual_texture/key.h"

#include <algorithm>
#include <atomic>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <limits>
#include <tbb/blocked_range.h>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for.h>
#include <vector>

namespace ybi
{

namespace
{
struct FeedbackAccumulator
{
    std::vector<unsigned long long> keys;
    uint32_t sampled = 0;
    uint32_t overflow = 0;
};

struct CPUIntegratorState
{
    const LaunchParams *params = nullptr;
    FeedbackAccumulator *feedback = nullptr;
    unsigned int pixelX = 0;
    unsigned int pixelY = 0;

    const LaunchParams &Params() const
    {
        return *params;
    }

    ybi::render::integrator::UInt2 LaunchIndex() const
    {
        return {pixelX, pixelY};
    }

    bool SampleTexture2D(const LaunchParams::MaterialTextureRef &textureRef,
                         float uu,
                         float vv,
                         ybi::render::integrator::Vec4 &outSample) const
    {
        if (textureRef.textureObject == 0ull || textureRef.valid == 0 || textureRef.width <= 0 ||
            textureRef.height <= 0)
        {
            return false;
        }

        bool blackS = false;
        bool blackT = false;
        const float wrappedU =
            ybi::render::integrator::ApplyWrapMode(uu, textureRef.wrapS, blackS);
        const float wrappedV =
            ybi::render::integrator::ApplyWrapMode(vv, textureRef.wrapT, blackT);
        if (blackS || blackT)
        {
            outSample = {0.0f, 0.0f, 0.0f, 0.0f};
            return true;
        }

        const int texelX = ybi::texture::TexelFromUnitUV(wrappedU, textureRef.width);
        const int texelY = ybi::texture::TexelFromUnitUV(wrappedV, textureRef.height);
        const int safeX =
            texelX < 0 ? 0 : (texelX >= textureRef.width ? textureRef.width - 1 : texelX);
        const int safeY =
            texelY < 0 ? 0 : (texelY >= textureRef.height ? textureRef.height - 1 : texelY);

        const unsigned char *pixels =
            reinterpret_cast<const unsigned char *>(textureRef.textureObject);
        const size_t index = (static_cast<size_t>(safeY) * static_cast<size_t>(textureRef.width) +
                              static_cast<size_t>(safeX)) *
                             4u;
        outSample = {pixels[index + 0] * (1.0f / 255.0f),
                     pixels[index + 1] * (1.0f / 255.0f),
                     pixels[index + 2] * (1.0f / 255.0f),
                     pixels[index + 3] * (1.0f / 255.0f)};
        return true;
    }

    void RecordFeedbackKey(unsigned long long key)
    {
        if (!feedback)
        {
            return;
        }

        feedback->sampled++;
        if (feedback->sampled <= static_cast<uint32_t>(params->feedbackCapacity))
        {
            feedback->keys.push_back(key);
        }
        else
        {
            feedback->overflow++;
        }
    }

    bool TraceOcclusion(const ybi::render::integrator::Vec3 &origin,
                        const ybi::render::integrator::Vec3 &direction,
                        float tMin,
                        float tMax) const
    {
        RTCRay shadowRay = {};
        shadowRay.org_x = origin.x;
        shadowRay.org_y = origin.y;
        shadowRay.org_z = origin.z;
        shadowRay.dir_x = direction.x;
        shadowRay.dir_y = direction.y;
        shadowRay.dir_z = direction.z;
        shadowRay.tnear = tMin;
        shadowRay.tfar = tMax;
        shadowRay.mask = ~0u;
        shadowRay.flags = 0u;

#if YBI_EMBREE_VERSION_MAJOR >= 4
        RTCRayQueryContext occlusionQueryContext;
        rtcInitRayQueryContext(&occlusionQueryContext);
        RTCOccludedArguments occludedArgs;
        rtcInitOccludedArguments(&occludedArgs);
        occludedArgs.context = &occlusionQueryContext;
        rtcOccluded1(reinterpret_cast<RTCScene>(params->traversable), &shadowRay, &occludedArgs);
#else
        RTCIntersectContext shadowContext;
        rtcInitIntersectContext(&shadowContext);
        rtcOccluded1(reinterpret_cast<RTCScene>(params->traversable), &shadowContext, &shadowRay);
#endif

        return shadowRay.tfar < 0.0f;
    }

    bool TraceLightOcclusion(const ybi::render::integrator::Vec3 &origin,
                             const ybi::render::integrator::Vec3 &direction,
                             float tMin,
                             float tMax,
                             int lightIndex) const;

    bool TraceClosest(const ybi::render::integrator::Vec3 &origin,
                      const ybi::render::integrator::Vec3 &direction,
                      float tMin,
                      float tMax,
                      ybi::render::integrator::HitInfo *outHit) const;

    void MaybeLogSampleSuccess() const
    {
        static std::atomic<unsigned int> logged{0};
        unsigned int expected = 0;
        if (logged.compare_exchange_strong(expected, 1u))
        {
            std::printf("TrySampleMaterialTexture succeeded\n");
        }
    }
};

static ybi::render::integrator::Vec3
TransformPoint3x4RowMajor(const float transform[12], const ybi::render::integrator::Vec3 &p)
{
    return ybi::render::integrator::Vec3(
        transform[0] * p.x + transform[1] * p.y + transform[2] * p.z + transform[3],
        transform[4] * p.x + transform[5] * p.y + transform[6] * p.z + transform[7],
        transform[8] * p.x + transform[9] * p.y + transform[10] * p.z + transform[11]);
}

static ybi::render::integrator::Vec3
TransformVector3x4RowMajor(const float transform[12], const ybi::render::integrator::Vec3 &v)
{
    return ybi::render::integrator::Vec3(
        transform[0] * v.x + transform[1] * v.y + transform[2] * v.z,
        transform[4] * v.x + transform[5] * v.y + transform[6] * v.z,
        transform[8] * v.x + transform[9] * v.y + transform[10] * v.z);
}

static ybi::render::integrator::Vec3
TransformNormal3x4RowMajor(const float transform[12], const ybi::render::integrator::Vec3 &n)
{
    const float m00 = transform[0];
    const float m01 = transform[1];
    const float m02 = transform[2];
    const float m10 = transform[4];
    const float m11 = transform[5];
    const float m12 = transform[6];
    const float m20 = transform[8];
    const float m21 = transform[9];
    const float m22 = transform[10];

    const float c00 = m11 * m22 - m12 * m21;
    const float c01 = m02 * m21 - m01 * m22;
    const float c02 = m01 * m12 - m02 * m11;
    const float c10 = m12 * m20 - m10 * m22;
    const float c11 = m00 * m22 - m02 * m20;
    const float c12 = m02 * m10 - m00 * m12;
    const float c20 = m10 * m21 - m11 * m20;
    const float c21 = m01 * m20 - m00 * m21;
    const float c22 = m00 * m11 - m01 * m10;
    const float det = m00 * c00 + m01 * c10 + m02 * c20;
    if (fabsf(det) <= 1.0e-12f)
    {
        return ybi::render::integrator::Normalize(TransformVector3x4RowMajor(transform, n));
    }
    const float invDet = 1.0f / det;
    const ybi::render::integrator::Vec3 transformed((c00 * n.x + c10 * n.y + c20 * n.z) * invDet,
                                                    (c01 * n.x + c11 * n.y + c21 * n.z) * invDet,
                                                    (c02 * n.x + c12 * n.y + c22 * n.z) * invDet);
    return ybi::render::integrator::Normalize(transformed);
}

static bool TryGetRayHitWorldTransform(const LaunchParams &params,
                                       const RTCRayHit &rayHit,
                                       float transform[12])
{
    unsigned int sceneGeomId = rayHit.hit.geomID;
    if (rayHit.hit.instID[0] != RTC_INVALID_GEOMETRY_ID)
    {
        sceneGeomId = rayHit.hit.instID[0];
    }

    rtcGetGeometryTransformFromScene(reinterpret_cast<RTCScene>(params.traversable),
                                     sceneGeomId,
                                     0.0f,
                                     RTC_FORMAT_FLOAT3X4_ROW_MAJOR,
                                     transform);
    return true;
}

static bool TryComputeTriangleWorldPositions(const LaunchParams &params,
                                             const RTCRayHit &rayHit,
                                             ybi::render::integrator::HitInfo *outHit)
{
    if (!outHit || params.instanceGeomRefs == 0ull || outHit->instanceId < 0 ||
        outHit->instanceId >= params.instanceGeomRefCount || outHit->primitiveIndex < 0)
    {
        return false;
    }

    ybi::render::integrator::Vec3 p0 = {};
    ybi::render::integrator::Vec3 p1 = {};
    ybi::render::integrator::Vec3 p2 = {};
    if (!ybi::TryGetTriangleLocalPositions(
            params, outHit->instanceId, outHit->primitiveIndex, p0, p1, p2))
    {
        return false;
    }

    float transform[12] = {};
    TryGetRayHitWorldTransform(params, rayHit, transform);

    outHit->worldTri0 = TransformPoint3x4RowMajor(transform, p0);
    outHit->worldTri1 = TransformPoint3x4RowMajor(transform, p1);
    outHit->worldTri2 = TransformPoint3x4RowMajor(transform, p2);
    outHit->hasWorldTriangle = true;
    return true;
}

static bool TracePrimary(const LaunchParams &params,
                         const ybi::render::integrator::Vec3 &origin,
                         const ybi::render::integrator::Vec3 &direction,
                         float tMin,
                         float tMax,
                         const ybi::render::integrator::RayDifferential *rayDiff,
                         ybi::render::integrator::HitInfo *outHit)
{
    RTCRayHit rayHit = {};
    rayHit.ray.org_x = origin.x;
    rayHit.ray.org_y = origin.y;
    rayHit.ray.org_z = origin.z;
    rayHit.ray.dir_x = direction.x;
    rayHit.ray.dir_y = direction.y;
    rayHit.ray.dir_z = direction.z;
    rayHit.ray.tnear = tMin;
    rayHit.ray.tfar = tMax;
    rayHit.ray.mask = ~0u;
    rayHit.ray.flags = 0u;
    rayHit.hit.geomID = RTC_INVALID_GEOMETRY_ID;

#if YBI_EMBREE_VERSION_MAJOR >= 4
    RTCRayQueryContext rayQueryContext;
    rtcInitRayQueryContext(&rayQueryContext);
    RTCIntersectArguments intersectArgs;
    rtcInitIntersectArguments(&intersectArgs);
    intersectArgs.context = &rayQueryContext;
    rtcIntersect1(reinterpret_cast<RTCScene>(params.traversable), &rayHit, &intersectArgs);
#else
    RTCIntersectContext context;
    rtcInitIntersectContext(&context);
    rtcIntersect1(reinterpret_cast<RTCScene>(params.traversable), &context, &rayHit);
#endif

    if (rayHit.hit.geomID == RTC_INVALID_GEOMETRY_ID)
    {
        return false;
    }

    if (outHit)
    {
        int instanceId = static_cast<int>(rayHit.hit.geomID);
        RTCScene embreeScene = reinterpret_cast<RTCScene>(params.traversable);
        const unsigned int instId = rayHit.hit.instID[0];
        const unsigned int geomId =
            (instId != RTC_INVALID_GEOMETRY_ID) ? instId : rayHit.hit.geomID;
        RTCGeometry geometry = rtcGetGeometry(embreeScene, geomId);
        if (geometry)
        {
            const void *userData = rtcGetGeometryUserData(geometry);
            if (userData)
            {
                const uint32_t refIndex =
                    static_cast<uint32_t>(reinterpret_cast<uintptr_t>(userData));
                if (refIndex != UINT32_MAX)
                {
                    instanceId = static_cast<int>(refIndex);
                }
            }
        }

        outHit->rayOrigin = origin;
        outHit->rayDir = direction;
        outHit->t = rayHit.ray.tfar;
        outHit->instanceId = instanceId;
        outHit->primitiveIndex = static_cast<int>(rayHit.hit.primID);

        float transform[12] = {};
        TryGetRayHitWorldTransform(params, rayHit, transform);
        outHit->geomNormal = TransformNormal3x4RowMajor(
            transform, Vec3(rayHit.hit.Ng_x, rayHit.hit.Ng_y, rayHit.hit.Ng_z));
        outHit->hasGeomNormal = true;

        const bool hasMeshRef =
            outHit->instanceId >= 0 && outHit->instanceId < params.instanceGeomRefCount;
        if (hasMeshRef)
        {
            const float u = rayHit.hit.u;
            const float v = rayHit.hit.v;
            outHit->barycentrics = ybi::render::integrator::Vec3(1.0f - u - v, u, v);
            outHit->hasBarycentrics = true;
            (void)TryComputeTriangleWorldPositions(params, rayHit, outHit);
            (void)ybi::TryComputeTriangleSurfacePartials(params, outHit);
            outHit->hasShadingNormal = ComputeTriangleShadingNormal(
                params, u, v, rayHit.hit.primID, outHit->instanceId, outHit->shadingNormal);

            if (outHit->hasShadingNormal)
            {
                outHit->shadingNormal =
                    TransformNormal3x4RowMajor(transform, outHit->shadingNormal);
            }
            if (rayDiff)
            {
                (void)ybi::TryComputeTriangleHitDifferentials(params, *rayDiff, outHit);
            }
        }
    }

    return true;
}

bool CPUIntegratorState::TraceClosest(const ybi::render::integrator::Vec3 &origin,
                                      const ybi::render::integrator::Vec3 &direction,
                                      float tMin,
                                      float tMax,
                                      ybi::render::integrator::HitInfo *outHit) const
{
    return TracePrimary(*params, origin, direction, tMin, tMax, nullptr, outHit);
}

bool CPUIntegratorState::TraceLightOcclusion(const ybi::render::integrator::Vec3 &origin,
                                             const ybi::render::integrator::Vec3 &direction,
                                             float tMin,
                                             float tMax,
                                             int lightIndex) const
{
    if (!ybi::render::integrator::LightHasShadowExcludes(*params, lightIndex))
    {
        return TraceOcclusion(origin, direction, tMin, tMax);
    }

    ybi::render::integrator::Vec3 currentOrigin = origin;
    float currentMax = tMax;
    constexpr int kMaxShadowSkips = 32;
    for (int skip = 0; skip < kMaxShadowSkips; ++skip)
    {
        ybi::render::integrator::HitInfo hit = {};
        if (!TracePrimary(*params, currentOrigin, direction, tMin, currentMax, nullptr, &hit))
        {
            return false;
        }
        if (!ybi::render::integrator::IsRefExcludedFromLightShadow(
                *params, lightIndex, hit.instanceId))
        {
            return true;
        }

        const float advance = hit.t + std::max(tMin, 1.0e-4f);
        if (advance >= currentMax)
        {
            return false;
        }
        currentOrigin = currentOrigin + direction * advance;
        currentMax -= advance;
    }

    return true;
}
} // namespace

bool CPUDispatchKernel(const DispatchParams &dispatchParams, RenderKernelId kernelId)
{
    if (dispatchParams.launchParamsDevice == 0 ||
        dispatchParams.launchParamsSize != sizeof(LaunchParams))
    {
        std::fprintf(stderr, "CPU kernel dispatch: launch params missing or invalid.\n");
        return false;
    }

    const LaunchParams *params =
        reinterpret_cast<const LaunchParams *>(dispatchParams.launchParamsDevice);
    if (!params || params->traversable == 0)
    {
        std::fprintf(stderr, "CPU kernel dispatch: traversable missing.\n");
        return false;
    }

    const uint32_t width = dispatchParams.width;
    const uint32_t height = dispatchParams.height;
    if (width == 0 || height == 0)
    {
        std::fprintf(stderr, "CPU kernel dispatch: invalid launch dims.\n");
        return false;
    }

    if (!dispatchParams.outputRGBA8.data() ||
        dispatchParams.outputRGBA8.numBytes() <
            static_cast<size_t>(width) * static_cast<size_t>(height) * 4u)
    {
        std::fprintf(stderr, "CPU kernel dispatch: output buffer missing.\n");
        return false;
    }

    tbb::enumerable_thread_specific<FeedbackAccumulator> feedbackTLS;

    tbb::parallel_for(
        tbb::blocked_range<uint32_t>(0, height), [&](const tbb::blocked_range<uint32_t> &range) {
            FeedbackAccumulator &local = feedbackTLS.local();
            CPUIntegratorState state = {};
            state.params = params;
            state.feedback = &local;

            for (uint32_t y = range.begin(); y != range.end(); ++y)
            {
                for (uint32_t x = 0; x < width; ++x)
                {
                    const size_t pixelIndex =
                        (static_cast<size_t>(y) * static_cast<size_t>(width) +
                         static_cast<size_t>(x)) *
                        4u;
                    uint8_t *dst = dispatchParams.outputRGBA8.data();
                    if (!ybi::render::integrator::ShouldRenderLaunchPixel(
                            params->singlePixelEnabled,
                            params->singlePixelX,
                            params->singlePixelY,
                            x,
                            y))
                    {
                        dst[pixelIndex + 0] = 0u;
                        dst[pixelIndex + 1] = 0u;
                        dst[pixelIndex + 2] = 0u;
                        dst[pixelIndex + 3] = 255u;
                        continue;
                    }

                    state.pixelX = x;
                    state.pixelY = y;

                    const ybi::render::integrator::RayDifferential rayDiff =
                        ybi::render::integrator::InitPerspectiveRayDifferential(
                            *params, static_cast<float>(x), static_cast<float>(y), width, height);
                    const ybi::render::integrator::Vec3 origin = rayDiff.origin;
                    const ybi::render::integrator::Vec3 direction = rayDiff.dir;

                    ybi::render::integrator::HitInfo hit = {};
                    const bool hitFound = TracePrimary(*params,
                                                       origin,
                                                       direction,
                                                       0.001f,
                                                       std::numeric_limits<float>::infinity(),
                                                       &rayDiff,
                                                       &hit);

                    unsigned int packed = 0u;
                    if (!hitFound)
                    {
                        const ybi::render::integrator::Vec3 sky =
                            ybi::render::integrator::SkyColor(direction);
                        packed = ybi::render::PackRGB8(sky.x, sky.y, sky.z);
                    }
                    else if (params->integrator == 2)
                    {
                        if (hit.hasBarycentrics)
                        {
                            ybi::render::integrator::IntegratorFeedbackOnly(state, hit);
                        }
                        packed = ybi::render::PackRGB8(0.0f, 0.0f, 0.0f);
                    }
                    else if (kernelId == RenderKernelId::AO)
                    {
                        packed = ybi::render::integrator::IntegratorAO(state, hit);
                    }
                    else if (kernelId == RenderKernelId::PathTrace)
                    {
                        packed =
                            ybi::render::integrator::IntegratorPathTrace(state, origin, direction);
                    }
                    else
                    {
                        packed = ybi::render::integrator::IntegratorPrimaryDiffuse(state, hit);
                    }

                    dst[pixelIndex + 0] = static_cast<uint8_t>(packed & 0xffu);
                    dst[pixelIndex + 1] = static_cast<uint8_t>((packed >> 8) & 0xffu);
                    dst[pixelIndex + 2] = static_cast<uint8_t>((packed >> 16) & 0xffu);
                    dst[pixelIndex + 3] = 255u;
                }
            }
        });

    if (params->feedbackStats != 0ull && params->feedbackKeys != 0ull &&
        params->feedbackCapacity > 0)
    {
        uint32_t totalSampled = 0;
        uint32_t totalOverflow = 0;
        for (FeedbackAccumulator &accum : feedbackTLS)
        {
            totalSampled += accum.sampled;
            totalOverflow += accum.overflow;
        }

        unsigned long long *keys = reinterpret_cast<unsigned long long *>(params->feedbackKeys);
        const size_t capacity = static_cast<size_t>(params->feedbackCapacity);
        size_t offset = 0;
        for (FeedbackAccumulator &accum : feedbackTLS)
        {
            if (offset >= capacity)
            {
                break;
            }
            const size_t copyCount =
                std::min(capacity - offset, static_cast<size_t>(accum.keys.size()));
            if (copyCount > 0)
            {
                std::memcpy(
                    keys + offset, accum.keys.data(), copyCount * sizeof(unsigned long long));
                offset += copyCount;
            }
        }

        unsigned int *stats = reinterpret_cast<unsigned int *>(params->feedbackStats);
        stats[0] = totalSampled;
        stats[1] = totalOverflow;
    }

    return true;
}

} // namespace ybi

#endif
