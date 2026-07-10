#pragma once 

#include "core/vec2f.h"
#include "core/vec3f.h"
#include "core/math.h"

namespace ybi
{

YBI_CPU_GPU uint32_t EncodeFloatToUint16(float f)
{
    return (uint32_t)((f + 1) / 2 * 65535.f + 0.5f);
}

YBI_CPU_GPU float DecodeUint16ToFloat(uint16_t u)
{
    return (float)(u / 65535.f * 2 - 1);
}

YBI_CPU_GPU float SignNotZero(float v)
{
    return v < 0.f ? -1.f : 1.f;
}

YBI_CPU_GPU uint32_t EncodeOctahedral(Vec3f v)
{
    Vec2f p(v.x, v.y);
    p = p * (1.f / (Abs(v.x) + Abs(v.y) + Abs(v.z)));
    if (v.z <= 0)
    {
        p = Vec2f((1.f - Abs(p.y)) * SignNotZero(p.x), (1.f - Abs(p.x)) * SignNotZero(p.y));
    }
    p.x = Clamp(p.x, -1.f, 1.f);
    p.y = Clamp(p.y, -1.f, 1.f);

    return (EncodeFloatToUint16(p[0]) << 16) | EncodeFloatToUint16(p[1]);
}

YBI_CPU_GPU Vec3f UnpackOctahedral(Vec2f v)
{
    Vec3f normal(v.x, v.y, 1 - Abs(v.x) - Abs(v.y));
    float t = Clamp(-normal.z, 0.f, 1.f);
    normal.x += normal.x >= 0.f ? -t : t;
    normal.y += normal.y >= 0.f ? -t : t;
    return Normalize(normal);
}

YBI_CPU_GPU Vec3f DecodeOctahedral(uint32_t packed)
{
    const uint16_t x = static_cast<uint16_t>(packed >> 16);
    const uint16_t y = static_cast<uint16_t>(packed & 0xffffu);
    return UnpackOctahedral(Vec2f(DecodeUint16ToFloat(x), DecodeUint16ToFloat(y)));
}

} // namespace ybi
