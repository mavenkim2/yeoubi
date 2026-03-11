#pragma once

#include "util/math_common.h"
#include <cmath>

namespace ybi
{

struct Vec3
{
    union
    {
        struct
        {
            float x;
            float y;
            float z;
        };
        struct
        {
            float r;
            float g;
            float b;
        };
    };

    Vec3() = default;
    YBI_INTEGRATOR_HD constexpr Vec3(float value) : x(value), y(value), z(value) {}
    YBI_INTEGRATOR_HD constexpr Vec3(float x_, float y_, float z_) : x(x_), y(y_), z(z_) {}

    YBI_INTEGRATOR_HD float operator[](int i) const { return *(&x + i); }
    YBI_INTEGRATOR_HD float &operator[](int i) { return *(&x + i); }
};

static_assert(sizeof(Vec3) == sizeof(float) * 3);

YBI_INTEGRATOR_HD Vec3 operator+(const Vec3 &a, const Vec3 &b)
{
    return Vec3(a.x + b.x, a.y + b.y, a.z + b.z);
}

YBI_INTEGRATOR_HD Vec3 &operator+=(Vec3 &a, const Vec3 &b)
{
    a.x += b.x;
    a.y += b.y;
    a.z += b.z;
    return a;
}

YBI_INTEGRATOR_HD Vec3 operator-(const Vec3 &a, const Vec3 &b)
{
    return Vec3(a.x - b.x, a.y - b.y, a.z - b.z);
}

YBI_INTEGRATOR_HD Vec3 &operator-=(Vec3 &a, const Vec3 &b)
{
    a.x -= b.x;
    a.y -= b.y;
    a.z -= b.z;
    return a;
}

YBI_INTEGRATOR_HD Vec3 operator-(const Vec3 &a)
{
    return Vec3(-a.x, -a.y, -a.z);
}

YBI_INTEGRATOR_HD Vec3 operator*(const Vec3 &a, const Vec3 &b)
{
    return Vec3(a.x * b.x, a.y * b.y, a.z * b.z);
}

YBI_INTEGRATOR_HD Vec3 &operator*=(Vec3 &a, const Vec3 &b)
{
    a.x *= b.x;
    a.y *= b.y;
    a.z *= b.z;
    return a;
}

YBI_INTEGRATOR_HD Vec3 operator*(const Vec3 &a, float b)
{
    return Vec3(a.x * b, a.y * b, a.z * b);
}

YBI_INTEGRATOR_HD Vec3 &operator*=(Vec3 &a, float b)
{
    a.x *= b;
    a.y *= b;
    a.z *= b;
    return a;
}

YBI_INTEGRATOR_HD Vec3 operator*(float a, const Vec3 &b)
{
    return b * a;
}

YBI_INTEGRATOR_HD Vec3 operator/(const Vec3 &a, const Vec3 &b)
{
    return Vec3(a.x / b.x, a.y / b.y, a.z / b.z);
}

YBI_INTEGRATOR_HD Vec3 operator/(const Vec3 &a, float b)
{
    return Vec3(a.x / b, a.y / b, a.z / b);
}

YBI_INTEGRATOR_HD Vec3 &operator/=(Vec3 &a, const Vec3 &b)
{
    a.x /= b.x;
    a.y /= b.y;
    a.z /= b.z;
    return a;
}

YBI_INTEGRATOR_HD Vec3 &operator/=(Vec3 &a, float b)
{
    a.x /= b;
    a.y /= b;
    a.z /= b;
    return a;
}

YBI_INTEGRATOR_HD float Dot(const Vec3 &a, const Vec3 &b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

YBI_INTEGRATOR_HD Vec3 Cross(const Vec3 &a, const Vec3 &b)
{
    return Vec3(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
}

YBI_INTEGRATOR_HD float LengthSquared(const Vec3 &a)
{
    return Dot(a, a);
}

YBI_INTEGRATOR_HD float Length(const Vec3 &a)
{
    return sqrtf(Dot(a, a));
}

YBI_INTEGRATOR_HD Vec3 Normalize(const Vec3 &a)
{
    const float lenSq = LengthSquared(a);
    if (lenSq <= 1e-20f)
    {
        return Vec3(0.0f, 0.0f, 1.0f);
    }
    return a * (1.0f / sqrtf(lenSq));
}

YBI_INTEGRATOR_HD Vec3 Clamp(const Vec3 &v, float lo, float hi)
{
    const float x = v.x < lo ? lo : (v.x > hi ? hi : v.x);
    const float y = v.y < lo ? lo : (v.y > hi ? hi : v.y);
    const float z = v.z < lo ? lo : (v.z > hi ? hi : v.z);
    return Vec3(x, y, z);
}

YBI_INTEGRATOR_HD Vec3 Lerp(const Vec3 &a, const Vec3 &b, float t)
{
    return a * (1.0f - t) + b * t;
}

YBI_INTEGRATOR_HD float MaxComponent(const Vec3 &v)
{
    return fmaxf(v.x, fmaxf(v.y, v.z));
}

YBI_INTEGRATOR_HD float Luminance(const Vec3 &v)
{
    return 0.2126f * v.x + 0.7152f * v.y + 0.0722f * v.z;
}

YBI_INTEGRATOR_HD Vec3 Reflect(const Vec3 &incident, const Vec3 &normal)
{
    return incident - 2.0f * Dot(incident, normal) * normal;
}

YBI_INTEGRATOR_HD Vec3 FaceForward(const Vec3 &normal, const Vec3 &referenceDirection)
{
    return Dot(normal, referenceDirection) < 0.0f ? normal : -normal;
}

YBI_INTEGRATOR_HD void BuildOrthonormalBasis(const Vec3 &n, Vec3 &t, Vec3 &b)
{
    const Vec3 up = fabsf(n.z) < 0.999f ? Vec3(0.0f, 0.0f, 1.0f) : Vec3(0.0f, 1.0f, 0.0f);
    t = Normalize(Cross(up, n));
    b = Normalize(Cross(n, t));
}

} // namespace ybi
