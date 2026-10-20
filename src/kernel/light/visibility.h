#pragma once

#include "kernel/light/types.h"
#include "kernel/raytracing/ray_offset.h"
#include "kernel/raytracing/trace.h"

namespace ybi
{

YBI_GPU bool IsLightSampleVisible(uint64_t bvhHandle,
                                  const LightSample &sample,
                                  const Vec3f &position,
                                  const Vec3f &geometricNormal)
{
    const Vec3f shadowNormal = Dot(sample.direction, geometricNormal) >= 0.f
                                  ? geometricNormal
                                  : -geometricNormal;
    const Vec3f shadowOrigin = OffsetRayOrigin(position, shadowNormal);
    Vec3f shadowDirection = sample.direction;
    float shadowTMax = 1.0e30f;
    if (sample.distance != FLOAT_INFINITY)
    {
        const Vec3f shadowToLight = sample.direction * sample.distance + (position - shadowOrigin);
        const float shadowDistance = Length(shadowToLight);
        if (shadowDistance == 0.f)
            return false;
        shadowDirection = shadowToLight / shadowDistance;
        shadowTMax = shadowDistance * (1.f - 1e-4f);
    }
    return TraceVisibility(bvhHandle, shadowOrigin, shadowDirection, shadowTMax);
}

} // namespace ybi
