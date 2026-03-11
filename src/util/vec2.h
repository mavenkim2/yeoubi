#pragma once

#include "util/math_common.h"
#include <cmath>

namespace ybi
{

struct Vec2
{
    union
    {
        struct
        {
            float x;
            float y;
        };
        struct
        {
            float r;
            float g;
        };
    };

    Vec2() = default;
    YBI_INTEGRATOR_HD constexpr Vec2(float value) : x(value), y(value) {}
    YBI_INTEGRATOR_HD constexpr Vec2(float x_, float y_) : x(x_), y(y_) {}

    YBI_INTEGRATOR_HD float operator[](int i) const { return *(&x + i); }
    YBI_INTEGRATOR_HD float &operator[](int i) { return *(&x + i); }
};

static_assert(sizeof(Vec2) == sizeof(float) * 2);

YBI_INTEGRATOR_HD Vec2 operator+(const Vec2 &a, const Vec2 &b)
{
    return Vec2(a.x + b.x, a.y + b.y);
}

YBI_INTEGRATOR_HD Vec2 &operator+=(Vec2 &a, const Vec2 &b)
{
    a.x += b.x;
    a.y += b.y;
    return a;
}

YBI_INTEGRATOR_HD Vec2 operator-(const Vec2 &a, const Vec2 &b)
{
    return Vec2(a.x - b.x, a.y - b.y);
}

YBI_INTEGRATOR_HD Vec2 operator-(const Vec2 &a)
{
    return Vec2(-a.x, -a.y);
}

YBI_INTEGRATOR_HD Vec2 &operator-=(Vec2 &a, const Vec2 &b)
{
    a.x -= b.x;
    a.y -= b.y;
    return a;
}

YBI_INTEGRATOR_HD Vec2 operator*(const Vec2 &a, const Vec2 &b)
{
    return Vec2(a.x * b.x, a.y * b.y);
}

YBI_INTEGRATOR_HD Vec2 &operator*=(Vec2 &a, const Vec2 &b)
{
    a.x *= b.x;
    a.y *= b.y;
    return a;
}

YBI_INTEGRATOR_HD Vec2 operator*(const Vec2 &a, float b)
{
    return Vec2(a.x * b, a.y * b);
}

YBI_INTEGRATOR_HD Vec2 &operator*=(Vec2 &a, float b)
{
    a.x *= b;
    a.y *= b;
    return a;
}

YBI_INTEGRATOR_HD Vec2 operator*(float a, const Vec2 &b)
{
    return b * a;
}

YBI_INTEGRATOR_HD Vec2 operator/(const Vec2 &a, const Vec2 &b)
{
    return Vec2(a.x / b.x, a.y / b.y);
}

YBI_INTEGRATOR_HD Vec2 operator/(const Vec2 &a, float b)
{
    return Vec2(a.x / b, a.y / b);
}

YBI_INTEGRATOR_HD Vec2 &operator/=(Vec2 &a, const Vec2 &b)
{
    a.x /= b.x;
    a.y /= b.y;
    return a;
}

YBI_INTEGRATOR_HD Vec2 &operator/=(Vec2 &a, float b)
{
    a.x /= b;
    a.y /= b;
    return a;
}

YBI_INTEGRATOR_HD Vec2 Lerp(const Vec2 &a, const Vec2 &b, float t)
{
    return a * (1.0f - t) + b * t;
}

YBI_INTEGRATOR_HD float Dot(const Vec2 &a, const Vec2 &b)
{
    return a.x * b.x + a.y * b.y;
}

YBI_INTEGRATOR_HD float Length(const Vec2 &a)
{
    return sqrtf(Dot(a, a));
}

} // namespace ybi
