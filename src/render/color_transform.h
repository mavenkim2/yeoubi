#pragma once

#include "util/math_common.h"
#include "util/vec3.h"

#include <cmath>

namespace ybi
{
namespace render
{

YBI_DEVICE float ToneMapPathRadianceComponent(float value)
{
    const float linear = value < 0.0f ? 0.0f : value;
    return 1.0f - expf(-linear);
}

YBI_DEVICE float LinearToSrgb(float value)
{
    if (value <= 0.0f)
    {
        return 0.0f;
    }
    if (value <= 0.0031308f)
    {
        return value * 12.92f;
    }
    return 1.055f * powf(value, 1.0f / 2.4f) - 0.055f;
}

YBI_DEVICE Vec3 DisplayMapPathRadiance(const Vec3 &radiance)
{
    return Vec3(LinearToSrgb(ToneMapPathRadianceComponent(radiance.x)),
                LinearToSrgb(ToneMapPathRadianceComponent(radiance.y)),
                LinearToSrgb(ToneMapPathRadianceComponent(radiance.z)));
}

} // namespace render
} // namespace ybi
