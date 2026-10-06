#pragma once

#include "core/lane.h"
#include <xmmintrin.h>

namespace ybi
{

struct Mask4
{
    __m128 value;
};

struct Float4
{
    static constexpr uint32_t Width = 4;
    __m128 value;

    Float4() = default;
    Float4(float scalar) : value(_mm_set1_ps(scalar)) {}
    explicit Float4(__m128 lanes) : value(lanes) {}

    void Store(float *values) const
    {
        _mm_storeu_ps(values, value);
    }

    Float4 &operator+=(Float4 other)
    {
        value = _mm_add_ps(value, other.value);
        return *this;
    }

    Float4 &operator/=(Float4 other)
    {
        value = _mm_div_ps(value, other.value);
        return *this;
    }
};

template <>
struct Vector3Traits<Float4>
{
    using Type = LaneVec3<Float4>;
};

inline Float4 LoadUnaligned(const float *values)
{
    return Float4(_mm_loadu_ps(values));
}

inline Vector3<Float4> LoadUnaligned(const Vec3f *values)
{
    static_assert(sizeof(Vec3f) == 3 * sizeof(float));
    // Four packed Vec3f entries occupy three SSE registers.
    const __m128 a = _mm_loadu_ps(&values[0].x); // x0 y0 z0 x1
    const __m128 b = _mm_loadu_ps(&values[1].y); // y1 z1 x2 y2
    const __m128 c = _mm_loadu_ps(&values[2].z); // z2 x3 y3 z3

    const __m128 x01 = _mm_shuffle_ps(a, a, _MM_SHUFFLE(3, 0, 3, 0));
    const __m128 x23 = _mm_shuffle_ps(b, c, _MM_SHUFFLE(1, 1, 2, 2));
    const __m128 y01 = _mm_shuffle_ps(a, b, _MM_SHUFFLE(0, 0, 1, 1));
    const __m128 y23 = _mm_shuffle_ps(b, c, _MM_SHUFFLE(2, 2, 3, 3));
    const __m128 z01 = _mm_shuffle_ps(a, b, _MM_SHUFFLE(1, 1, 2, 2));
    const __m128 z23 = _mm_shuffle_ps(c, c, _MM_SHUFFLE(3, 0, 3, 0));
    return {Float4(_mm_shuffle_ps(x01, x23, _MM_SHUFFLE(2, 0, 1, 0))),
            Float4(_mm_shuffle_ps(y01, y23, _MM_SHUFFLE(2, 0, 2, 0))),
            Float4(_mm_shuffle_ps(z01, z23, _MM_SHUFFLE(1, 0, 2, 0)))};
}

inline Float4 operator+(Float4 a, Float4 b)
{
    return Float4(_mm_add_ps(a.value, b.value));
}

inline Float4 operator-(Float4 a, Float4 b)
{
    return Float4(_mm_sub_ps(a.value, b.value));
}

inline Float4 operator*(Float4 a, Float4 b)
{
    return Float4(_mm_mul_ps(a.value, b.value));
}

inline Float4 operator/(Float4 a, Float4 b)
{
    return Float4(_mm_div_ps(a.value, b.value));
}

inline Mask4 operator==(Float4 a, Float4 b)
{
    return {_mm_cmpeq_ps(a.value, b.value)};
}

inline Mask4 operator>(Float4 a, Float4 b)
{
    return {_mm_cmpgt_ps(a.value, b.value)};
}

inline Mask4 operator>=(Float4 a, Float4 b)
{
    return {_mm_cmpge_ps(a.value, b.value)};
}

inline Float4 Select(Mask4 mask, Float4 ifTrue, Float4 ifFalse)
{
    return Float4(
        _mm_or_ps(_mm_and_ps(mask.value, ifTrue.value), _mm_andnot_ps(mask.value, ifFalse.value)));
}

inline bool All(Mask4 mask)
{
    return _mm_movemask_ps(mask.value) == 15;
}

inline float ReduceSum(Float4 lanes)
{
    const __m128 pairs = _mm_add_ps(lanes.value, _mm_movehl_ps(lanes.value, lanes.value));
    const __m128 other = _mm_shuffle_ps(pairs, pairs, _MM_SHUFFLE(1, 1, 1, 1));
    return _mm_cvtss_f32(_mm_add_ss(pairs, other));
}

// Initial implementation: arithmetic uses SIMD; exp and log use scalar libm per lane.
inline Float4 Exp(Float4 value)
{
    float lanes[Float4::Width];
    value.Store(lanes);
    for (float &lane : lanes)
    {
        lane = expf(lane);
    }
    return LoadUnaligned(lanes);
}

inline Float4 Log(Float4 value)
{
    float lanes[Float4::Width];
    value.Store(lanes);
    for (float &lane : lanes)
    {
        lane = logf(lane);
    }
    return LoadUnaligned(lanes);
}

} // namespace ybi
