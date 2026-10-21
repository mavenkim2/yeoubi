#pragma once 

#include "core/vec2f.h"
#include "core/vec3f.h"
#include "core/math.h"
#include "core/simd_float4.h"

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

#if defined(SIMD_4)
// Decode four packed directions into coordinate lanes. Reads exactly four uint32_t values.
inline Vector3<Float4> DecodeOctahedral4(const uint32_t *values)
{
    const __m128i packed = _mm_loadu_si128(reinterpret_cast<const __m128i *>(values));
    const __m128 zero = _mm_setzero_ps();
    const __m128 one = _mm_set1_ps(1.f);
    const __m128 two = _mm_set1_ps(2.f);
    const __m128 scale = _mm_set1_ps(65535.f);
    const __m128 sign = _mm_set1_ps(-0.f);
    __m128 x = _mm_cvtepi32_ps(_mm_srli_epi32(packed, 16));
    __m128 y = _mm_cvtepi32_ps(_mm_and_si128(packed, _mm_set1_epi32(0xffff)));
    x = _mm_sub_ps(_mm_mul_ps(_mm_div_ps(x, scale), two), one);
    y = _mm_sub_ps(_mm_mul_ps(_mm_div_ps(y, scale), two), one);
    const __m128 z = _mm_sub_ps(_mm_sub_ps(one, _mm_andnot_ps(sign, x)),
                               _mm_andnot_ps(sign, y));
    const __m128 t = _mm_min_ps(_mm_max_ps(_mm_sub_ps(zero, z), zero), one);
    const __m128 negativeT = _mm_sub_ps(zero, t);
    x = _mm_add_ps(x, _mm_blendv_ps(t, negativeT, _mm_cmpge_ps(x, zero)));
    y = _mm_add_ps(y, _mm_blendv_ps(t, negativeT, _mm_cmpge_ps(y, zero)));
    const __m128 lengthSquared = _mm_add_ps(_mm_add_ps(_mm_mul_ps(x, x), _mm_mul_ps(y, y)),
                                           _mm_mul_ps(z, z));
    const __m128 invLength = _mm_div_ps(one, _mm_sqrt_ps(lengthSquared));
    return {Float4(_mm_mul_ps(x, invLength)),
            Float4(_mm_mul_ps(y, invLength)),
            Float4(_mm_mul_ps(z, invLength))};
}
#endif

} // namespace ybi
