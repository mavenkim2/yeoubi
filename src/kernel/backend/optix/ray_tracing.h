#pragma once

#include "core/vec3f.h"
#include "device/cuda_compat.h"
#include "kernel/backend/optix/visibility_mask.h"
#include "kernel/state/types.h"
#include <optix_device.h>

namespace ybi
{

YBI_GPU IntersectionInfo TraceOptixRay(uint64_t bvhHandle,
                                      const Vec3f &rayOrigin,
                                      const Vec3f &rayDirection,
                                      float tMin,
                                      float tMax)
{
    IntersectionInfo info{};
    const uint64_t ptr = reinterpret_cast<uint64_t>(&info);
    uint32_t low = static_cast<uint32_t>(ptr);
    uint32_t high = static_cast<uint32_t>(ptr >> 32u);

    optixTrace(static_cast<OptixTraversableHandle>(bvhHandle),
               make_float3(rayOrigin.x, rayOrigin.y, rayOrigin.z),
               make_float3(rayDirection.x, rayDirection.y, rayDirection.z),
               tMin,
               tMax,
               0.f,
               OPTIX_VISIBILITY_ALL,
               OPTIX_RAY_FLAG_NONE,
               0,
               0,
               0,
               low,
               high);

    return info;
}

YBI_GPU bool TraceOptixVisibility(uint64_t bvhHandle,
                                  const Vec3f &rayOrigin,
                                  const Vec3f &rayDirection,
                                  float tMax)
{
    if (tMax <= 0.f)
        return false;

    IntersectionInfo shadowInfo{};
    const uint64_t ptr = reinterpret_cast<uint64_t>(&shadowInfo);
    uint32_t low = static_cast<uint32_t>(ptr);
    uint32_t high = static_cast<uint32_t>(ptr >> 32u);

    optixTrace(static_cast<OptixTraversableHandle>(bvhHandle),
               make_float3(rayOrigin.x, rayOrigin.y, rayOrigin.z),
               make_float3(rayDirection.x, rayDirection.y, rayDirection.z),
               0.f,
               tMax,
               0.f,
               OPTIX_VISIBILITY_SHADOW,
               OPTIX_RAY_FLAG_TERMINATE_ON_FIRST_HIT,
               0,
               0,
               0,
               low,
               high);

    return shadowInfo.type == INTERSECTION_TYPE_NONE;
}

} // namespace ybi
