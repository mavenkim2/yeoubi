#pragma once

#include "util/assert.h"
#include "util/basemath.h"
#include "util/forceinline.h"
#include <cmath>

YBI_NAMESPACE_BEGIN

#ifndef GPU_INTEGRATOR
struct float3
{
    float x, y, z;
    __forceinline float operator[](int i) const
    {
        YBI_ASSERT(i >= 0);
        YBI_ASSERT(i < 3);
        return *(&x + i);
    }
    __forceinline float &operator[](int i)
    {
        YBI_ASSERT(i >= 0);
        YBI_ASSERT(i < 3);
        return *(&x + i);
    }
};

__forceinline float3 make_float3(const float x, const float y, const float z)
{
    return {x, y, z};
}
#endif

#ifndef YBI_INTEGRATOR_HD
#define YBI_INTEGRATOR_HD inline
#endif

__forceinline float3 make_float3(const float a)
{
    return {a, a, a};
}


// operators ////////////////////////////////////////////////////////////////

__forceinline float3 operator+(const float3 &a, const float3 &b)
{
    return {a.x + b.x, a.y + b.y, a.z + b.z};
}

__forceinline float3 &operator-=(float3 &a, const float3 &b)
{
    a.x -= b.x;
    a.y -= b.y;
    a.z -= b.z;
    return a;
}

__forceinline float3 &operator+=(float3 &a, const float3 &b)
{
    a.x += b.x;
    a.y += b.y;
    a.z += b.z;
    return a;
}

__forceinline float3 operator-(const float3 &a, const float3 &b)
{
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}

__forceinline float3 operator*(const float3 &a, const float3 &b)
{
    return {a.x * b.x, a.y * b.y, a.z * b.z};
}

__forceinline float3 operator*(const float3 &a, const float b)
{
    return {a.x * b, a.y * b, a.z * b};
}

__forceinline float3 operator*(float a, const float3 &b)
{
    return b * a;
}

__forceinline float3 operator/(const float3 &a, const float3 &b)
{
    return {a.x / b.x, a.y / b.y, a.z / b.z};
}

__forceinline float3 operator/(const float3 &a, const float b)
{
    return {a.x / b, a.y / b, a.z / b};
}

__forceinline float dot(const float3 &a, const float3 &b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

__forceinline float3 cross(const float3 &a, const float3 &b)
{
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}

__forceinline float length_squared(const float3 &a)
{
    return dot(a, a);
}

__forceinline float length(const float3 &a)
{
    return std::sqrt(dot(a, a));
}

__forceinline float3 normalize(const float3 &a)
{
    const float lenSq = length_squared(a);
    if (lenSq <= 1e-20f)
    {
        return make_float3(0.0f, 0.0f, 1.0f);
    }
    const float invLen = 1.0f / std::sqrt(lenSq);
    return a * invLen;
}

YBI_NAMESPACE_END
