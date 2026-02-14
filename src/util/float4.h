#pragma once

#include "util/assert.h"

YBI_NAMESPACE_BEGIN

struct float4
{
    float x, y, z, w;
    __forceinline float operator[](int i) const
    {
        YBI_ASSERT(i >= 0);
        YBI_ASSERT(i < 4);
        return *(&x + i);
    }
    __forceinline float &operator[](int i)
    {
        YBI_ASSERT(i >= 0);
        YBI_ASSERT(i < 4);
        return *(&x + i);
    }
};

__forceinline float4 make_float4(const float x, const float y, const float z, const float w)
{
    return {x, y, z, w};
}

// operator overloads ////////////////////////////////////////////////////////

__forceinline float4 operator+(const float4 &a, const float4 &b)
{
    return {a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w};
}

__forceinline float4 operator-(const float4 &a, const float4 &b)
{
    return {a.x - b.x, a.y - b.y, a.z - b.z, a.w - b.w};
}

__forceinline float4 operator*(const float4 &a, const float4 &b)
{
    return {a.x * b.x, a.y * b.y, a.z * b.z, a.w * b.w};
}

__forceinline float4 operator*(const float4 &a, const float b)
{
    return {a.x * b, a.y * b, a.z * b, a.w * b};
}

__forceinline float4 operator*(const float b, const float4 &a)
{
    return a * b;
}
// utility functions //////////////////////////////////////////////////////////

__forceinline float4 lerp(const float4 &a, const float4 &b, const float t)
{
    return a * (1.f - t) + b * t;
}  
__forceinline float dot(const float4 &a, const float4 &b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
}
__forceinline float length(const float4 &a)
{
    return sqrtf(dot(a, a));
}

YBI_NAMESPACE_END
