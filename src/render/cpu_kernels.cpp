#include "render/cpu_kernels.h"

#if defined(WITH_EMBREE)

#include "device/cpu_device.h"
#include "render/integrator_ao.h"
#include "render/integrator_primary.h"
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

YBI_NAMESPACE_BEGIN

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
        const int safeX = texelX < 0 ? 0 : (texelX >= textureRef.width ? textureRef.width - 1 : texelX);
        const int safeY = texelY < 0 ? 0 : (texelY >= textureRef.height ? textureRef.height - 1 : texelY);

        const unsigned char *pixels =
            reinterpret_cast<const unsigned char *>(textureRef.textureObject);
        const size_t index =
            (static_cast<size_t>(safeY) * static_cast<size_t>(textureRef.width) +
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

static ybi::render::integrator::Vec3 ComputeDirection(const LaunchParams &params,
                                                      unsigned int x,
                                                      unsigned int y,
                                                      unsigned int width,
                                                      unsigned int height)
{
    const float ndcX = (float(x) + 0.5f) / float(width) * 2.0f - 1.0f;
    const float ndcY = 1.0f - (float(y) + 0.5f) / float(height) * 2.0f;
    const ybi::render::integrator::Vec3 dir = ybi::render::integrator::Normalize(
        ybi::render::integrator::MakeVec3(params.cameraU.x * ndcX + params.cameraV.x * ndcY +
                                              params.cameraW.x,
                                          params.cameraU.y * ndcX + params.cameraV.y * ndcY +
                                              params.cameraW.y,
                                          params.cameraU.z * ndcX + params.cameraV.z * ndcY +
                                              params.cameraW.z));
    return dir;
}

static bool TracePrimary(const LaunchParams &params,
                         const ybi::render::integrator::Vec3 &origin,
                         const ybi::render::integrator::Vec3 &direction,
                         ybi::render::integrator::HitInfo *outHit)
{
    RTCRayHit rayHit = {};
    rayHit.ray.org_x = origin.x;
    rayHit.ray.org_y = origin.y;
    rayHit.ray.org_z = origin.z;
    rayHit.ray.dir_x = direction.x;
    rayHit.ray.dir_y = direction.y;
    rayHit.ray.dir_z = direction.z;
    rayHit.ray.tnear = 0.001f;
    rayHit.ray.tfar = std::numeric_limits<float>::infinity();
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
        outHit->geomNormal = ybi::render::integrator::Normalize(
            ybi::render::integrator::MakeVec3(rayHit.hit.Ng_x, rayHit.hit.Ng_y, rayHit.hit.Ng_z));
        outHit->hasGeomNormal = true;

        const bool hasMeshRef =
            outHit->instanceId >= 0 && outHit->instanceId < params.instanceGeomRefCount;
        if (hasMeshRef)
        {
            const float u = rayHit.hit.u;
            const float v = rayHit.hit.v;
            outHit->barycentrics = ybi::render::integrator::MakeVec3(1.0f - u - v, u, v);
            outHit->hasBarycentrics = true;
        }
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

    tbb::parallel_for(tbb::blocked_range<uint32_t>(0, height),
                      [&](const tbb::blocked_range<uint32_t> &range) {
                          FeedbackAccumulator &local = feedbackTLS.local();
                          CPUIntegratorState state = {};
                          state.params = params;
                          state.feedback = &local;

                          for (uint32_t y = range.begin(); y != range.end(); ++y)
                          {
                              for (uint32_t x = 0; x < width; ++x)
                              {
                                  state.pixelX = x;
                                  state.pixelY = y;

                                  const ybi::render::integrator::Vec3 origin =
                                      ybi::render::integrator::MakeVec3(params->cameraOrigin.x,
                                                                        params->cameraOrigin.y,
                                                                        params->cameraOrigin.z);
                                  const ybi::render::integrator::Vec3 direction =
                                      ComputeDirection(*params, x, y, width, height);

                                  ybi::render::integrator::HitInfo hit = {};
                                  const bool hitFound = TracePrimary(*params, origin, direction, &hit);

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
                                  else
                                  {
                                      packed =
                                          ybi::render::integrator::IntegratorPrimaryDiffuse(state, hit);
                                  }

                                  const size_t pixelIndex =
                                      (static_cast<size_t>(y) * static_cast<size_t>(width) +
                                       static_cast<size_t>(x)) *
                                      4u;
                                  uint8_t *dst = dispatchParams.outputRGBA8.data();
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

        unsigned long long *keys =
            reinterpret_cast<unsigned long long *>(params->feedbackKeys);
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
                std::memcpy(keys + offset, accum.keys.data(), copyCount * sizeof(unsigned long long));
                offset += copyCount;
            }
        }

        unsigned int *stats = reinterpret_cast<unsigned int *>(params->feedbackStats);
        stats[0] = totalSampled;
        stats[1] = totalOverflow;
    }

    return true;
}

YBI_NAMESPACE_END

#endif
