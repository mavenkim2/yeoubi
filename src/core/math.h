#pragma once

#include "device/cuda_compat.h"
#include <cmath>

#if defined(YBI_HAS_SSE)
#include <xmmintrin.h>
#endif

#define YBI_PI 3.14159265358979323846f
#define YBI_RCP_PI 0.318309886f
#define YBI_PI_OVER_2 1.57079632679489661923f
#define YBI_PI_OVER_4 0.78539816339744830962f
#define FLOAT_BEFORE_ONE 0x1.fffffep-1f

namespace ybi
{

template <typename T>
YBI_CPU_GPU constexpr T Clamp(T value, T minValue, T maxValue)
{
    return value < minValue ? minValue : (value > maxValue ? maxValue : value);
}

template <typename T>
YBI_CPU_GPU constexpr T Min(T a, T b)
{
    return a < b ? a : b;
}

template <typename T>
YBI_CPU_GPU constexpr T Max(T a, T b)
{
    return a > b ? a : b;
}

YBI_CPU_GPU float Floor(float a)
{
    return floorf(a);
}

YBI_CPU_GPU float Abs(float a)
{
    return a < 0.f ? -a : a;
}

YBI_CPU_GPU float Cos(float a)
{
    return cosf(a);
}

YBI_CPU_GPU float Sin(float a)
{
    return sinf(a);
}

// Returns an approximate reciprocal square root when the target provides a
// hardware implementation. Other targets use the precise scalar expression.
YBI_CPU_GPU float RcpSqrt(float a)
{
#if defined(__CUDA_ARCH__)
    return rsqrtf(a);
#elif defined(YBI_HAS_SSE)
    return _mm_cvtss_f32(_mm_rsqrt_ss(_mm_set_ss(a)));
#else
    return 1.0f / sqrtf(a);
#endif
}

YBI_CPU_GPU constexpr float DegreesToRadians(float degrees)
{
    return degrees * (YBI_PI / 180.f);
}

template <typename T>
YBI_CPU_GPU void Swap(T &a, T &b)
{
    T temp = a;
    a = b;
    b = temp;
}

template <typename T>
YBI_CPU_GPU T Sqr(const T &v)
{
    return v * v;
}

} // namespace ybi
