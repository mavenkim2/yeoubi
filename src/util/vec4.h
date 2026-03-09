#pragma once

#include "util/math_common.h"
#include <cmath>

YBI_NAMESPACE_BEGIN

struct Vec4
{
    float x;
    float y;
    float z;
    float w;

    Vec4() = default;
    YBI_INTEGRATOR_HD constexpr Vec4(float value) : x(value), y(value), z(value), w(value) {}
    YBI_INTEGRATOR_HD constexpr Vec4(float x_, float y_, float z_, float w_)
        : x(x_), y(y_), z(z_), w(w_)
    {
    }

    YBI_INTEGRATOR_HD float operator[](int i) const { return *(&x + i); }
    YBI_INTEGRATOR_HD float &operator[](int i) { return *(&x + i); }
};

static_assert(sizeof(Vec4) == sizeof(float) * 4);

YBI_INTEGRATOR_HD Vec4 operator+(const Vec4 &a, const Vec4 &b)
{
    return Vec4(a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w);
}

YBI_INTEGRATOR_HD Vec4 &operator+=(Vec4 &a, const Vec4 &b)
{
    a.x += b.x;
    a.y += b.y;
    a.z += b.z;
    a.w += b.w;
    return a;
}

YBI_INTEGRATOR_HD Vec4 operator-(const Vec4 &a, const Vec4 &b)
{
    return Vec4(a.x - b.x, a.y - b.y, a.z - b.z, a.w - b.w);
}

YBI_INTEGRATOR_HD Vec4 operator*(const Vec4 &a, const Vec4 &b)
{
    return Vec4(a.x * b.x, a.y * b.y, a.z * b.z, a.w * b.w);
}

YBI_INTEGRATOR_HD Vec4 operator*(const Vec4 &a, float b)
{
    return Vec4(a.x * b, a.y * b, a.z * b, a.w * b);
}

YBI_INTEGRATOR_HD Vec4 operator*(float b, const Vec4 &a)
{
    return a * b;
}

YBI_INTEGRATOR_HD Vec4 Lerp(const Vec4 &a, const Vec4 &b, float t)
{
    return a * (1.0f - t) + b * t;
}

YBI_INTEGRATOR_HD float Dot(const Vec4 &a, const Vec4 &b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
}

YBI_INTEGRATOR_HD float Length(const Vec4 &a)
{
    return sqrtf(Dot(a, a));
}

YBI_NAMESPACE_END
