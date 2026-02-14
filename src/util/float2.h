#pragma once

#include "util/assert.h"

YBI_NAMESPACE_BEGIN

struct float2
{
    float x, y;
    __forceinline float operator[](int i) const
    {
        YBI_ASSERT(i >= 0);
        YBI_ASSERT(i < 2);
        return *(&x + i);
    }
    __forceinline float &operator[](int i)
    {
        YBI_ASSERT(i >= 0);
        YBI_ASSERT(i < 2);
        return *(&x + i);
    }
};

__forceinline float2 make_float2(const float x, const float y)
{
    return {x, y};
}

__forceinline float2 make_float2(const float a)
{
    return {a, a};
}

// operators ////////////////////////////////////////////////////////////////

__forceinline float2 operator+(const float2 &a, const float2 &b)
{
    return {a.x + b.x, a.y + b.y};
}

__forceinline float2 &operator+=(float2 &a, const float2 &b)
{
    a.x += b.x;
    a.y += b.y;
    return a;
}

__forceinline float2 operator-(const float2 &a, const float2 &b)
{
    return {a.x - b.x, a.y - b.y};
}

__forceinline float2 operator*(const float2 &a, const float2 &b)
{
    return {a.x * b.x, a.y * b.y};
}

__forceinline float2 operator/(const float2 &a, const float2 &b)
{
    return {a.x / b.x, a.y / b.y};
}

__forceinline float2 operator/(const float2 &a, const float b)
{
    return {a.x / b, a.y / b};
}

__forceinline float2 operator*(const float2 &a, const float b)
{
    return {a.x * b, a.y * b};
}

// utility functions ////////////////////////////////////////////////////////

__forceinline float2 lerp(const float2 &a, const float2 &b, const float t)
{
    return a * (1.f - t) + b * t;
}

__forceinline float dot(const float2 &a, const float2 &b)
{
    return a.x * b.x + a.y * b.y;
}

__forceinline float length(const float2 &a)
{
    return sqrt(dot(a, a));
}
YBI_NAMESPACE_END
