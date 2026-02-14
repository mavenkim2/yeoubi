#pragma once

#include "util/assert.h"
#include <cmath>

YBI_NAMESPACE_BEGIN

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

__forceinline float3 make_float3(const float a)
{
    return {a, a, a};
}

// operators ////////////////////////////////////////////////////////////////

__forceinline float3 operator+(const float3 &a, const float3 &b)
{
    return {a.x + b.x, a.y + b.y, a.z + b.z};
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

__forceinline float length(const float3 &a)
{
    return sqrt(dot(a, a));
}

YBI_NAMESPACE_END
