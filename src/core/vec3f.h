#pragma once

#include "core/assert.h"
#include "core/math.h"
#include "core/vec2f.h"
#include "device/cuda_compat.h"
#include <cmath>
#include <cstddef>

namespace ybi
{

class Vec3f
{
public:
    float x;
    float y;
    float z;

    Vec3f() = default;
    YBI_CPU_GPU constexpr explicit Vec3f(float value) : x(value), y(value), z(value) {}
    YBI_CPU_GPU constexpr explicit Vec3f(const Vec2f &v, float value) : x(v.x), y(v.y), z(value) {}
    YBI_CPU_GPU constexpr explicit Vec3f(float value, const Vec2f &v) : x(value), y(v.x), z(v.y) {}
    YBI_CPU_GPU constexpr Vec3f(float x, float y, float z) : x(x), y(y), z(z) {}

    YBI_CPU_GPU float *Data()
    {
        return &x;
    }

    YBI_CPU_GPU const float *Data() const
    {
        return &x;
    }

    YBI_CPU_GPU float &operator[](size_t index)
    {
        YBI_ASSERT(index < 3);
        return Data()[index];
    }

    YBI_CPU_GPU const float &operator[](size_t index) const
    {
        YBI_ASSERT(index < 3);
        return Data()[index];
    }

    YBI_CPU_GPU constexpr Vec2f Xy() const
    {
        return Vec2f(x, y);
    }

    YBI_CPU_GPU constexpr Vec3f operator+() const
    {
        return *this;
    }

    YBI_CPU_GPU constexpr Vec3f operator-() const
    {
        return Vec3f(-x, -y, -z);
    }

    YBI_CPU_GPU Vec3f &operator+=(Vec3f rhs)
    {
        x += rhs.x;
        y += rhs.y;
        z += rhs.z;
        return *this;
    }

    YBI_CPU_GPU Vec3f &operator-=(Vec3f rhs)
    {
        x -= rhs.x;
        y -= rhs.y;
        z -= rhs.z;
        return *this;
    }

    YBI_CPU_GPU Vec3f &operator*=(Vec3f rhs)
    {
        x *= rhs.x;
        y *= rhs.y;
        z *= rhs.z;
        return *this;
    }

    YBI_CPU_GPU Vec3f &operator*=(float scalar)
    {
        x *= scalar;
        y *= scalar;
        z *= scalar;
        return *this;
    }

    YBI_CPU_GPU Vec3f &operator/=(float scalar)
    {
        YBI_ASSERT(scalar != 0.0f);
        const float invScalar = 1.0f / scalar;
        return *this *= invScalar;
    }

    YBI_CPU_GPU Vec3f &operator/=(Vec3f rhs)
    {
        YBI_ASSERT(rhs.x != 0.0f);
        YBI_ASSERT(rhs.y != 0.0f);
        YBI_ASSERT(rhs.z != 0.0f);
        x /= rhs.x;
        y /= rhs.y;
        z /= rhs.z;
        return *this;
    }

    YBI_CPU_GPU float MaxComponent() const
    {
        return Max<float>(x, Max<float>(y, z));
    }

    YBI_CPU_GPU float MinComponent() const
    {
        return Min<float>(x, Min<float>(y, z));
    }
};

YBI_CPU_GPU constexpr Vec3f operator+(Vec3f lhs, Vec3f rhs)
{
    return Vec3f(lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z);
}

YBI_CPU_GPU constexpr Vec3f operator+(Vec3f lhs, float rhs)
{
    return Vec3f(lhs.x + rhs, lhs.y + rhs, lhs.z + rhs);
}

YBI_CPU_GPU constexpr Vec3f operator+(float lhs, Vec3f rhs)
{
    return rhs + lhs;
}

YBI_CPU_GPU constexpr Vec3f operator-(Vec3f lhs, Vec3f rhs)
{
    return Vec3f(lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z);
}

YBI_CPU_GPU constexpr Vec3f operator-(Vec3f lhs, float rhs)
{
    return Vec3f(lhs.x - rhs, lhs.y - rhs, lhs.z - rhs);
}

YBI_CPU_GPU constexpr Vec3f operator-(float lhs, Vec3f rhs)
{
    return Vec3f(lhs - rhs.x, lhs - rhs.y, lhs - rhs.z);
}

YBI_CPU_GPU constexpr Vec3f operator*(Vec3f lhs, Vec3f rhs)
{
    return Vec3f(lhs.x * rhs.x, lhs.y * rhs.y, lhs.z * rhs.z);
}

YBI_CPU_GPU constexpr Vec3f operator*(Vec3f vector, float scalar)
{
    return Vec3f(vector.x * scalar, vector.y * scalar, vector.z * scalar);
}

YBI_CPU_GPU constexpr Vec3f operator*(float scalar, Vec3f vector)
{
    return vector * scalar;
}

YBI_CPU_GPU Vec3f operator/(Vec3f vector, float scalar)
{
    YBI_ASSERT(scalar != 0.0f);
    const float invScalar = 1.0f / scalar;
    return vector * invScalar;
}

YBI_CPU_GPU Vec3f operator/(float scalar, Vec3f v)
{
    return Vec3f(scalar / v.x, scalar / v.y, scalar / v.z);
}

YBI_CPU_GPU Vec3f operator/(Vec3f lhs, Vec3f rhs)
{
    lhs /= rhs;
    return lhs;
}

YBI_CPU_GPU constexpr bool operator==(Vec3f lhs, Vec3f rhs)
{
    return lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z;
}

YBI_CPU_GPU constexpr bool operator!=(Vec3f lhs, Vec3f rhs)
{
    return !(lhs == rhs);
}

YBI_CPU_GPU constexpr float Dot(Vec3f lhs, Vec3f rhs)
{
    return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z;
}

YBI_CPU_GPU constexpr Vec3f Cross(Vec3f lhs, Vec3f rhs)
{
    return Vec3f(lhs.y * rhs.z - lhs.z * rhs.y,
                 lhs.z * rhs.x - lhs.x * rhs.z,
                 lhs.x * rhs.y - lhs.y * rhs.x);
}

YBI_CPU_GPU constexpr float LengthSquared(Vec3f vector)
{
    return Dot(vector, vector);
}

YBI_CPU_GPU float Length(Vec3f vector)
{
    return std::sqrt(LengthSquared(vector));
}

YBI_CPU_GPU Vec3f Normalize(Vec3f vector)
{
    const float length = Length(vector);
    YBI_ERROR(length != 0.0f, "%f %f %f\n", vector.x, vector.y, vector.z);
    return vector / length;
}

YBI_CPU_GPU Vec3f Sqrt(Vec3f v)
{
    return Vec3f(Sqrt(v.x), Sqrt(v.y), Sqrt(v.z));
}

YBI_CPU_GPU constexpr Vec3f Lerp(Vec3f a, Vec3f b, float t)
{
    return a + (b - a) * t;
}

YBI_CPU_GPU constexpr Vec3f Lerp(float a, Vec3f b, float t)
{
    return a + (b - a) * t;
}

YBI_CPU_GPU constexpr Vec3f Lerp(Vec3f a, float b, float t)
{
    return a + (b - a) * t;
}

YBI_CPU_GPU constexpr Vec3f Min(Vec3f a, Vec3f b)
{
    return Vec3f(Min<float>(a.x, b.x), Min<float>(a.y, b.y), Min<float>(a.z, b.z));
}

YBI_CPU_GPU constexpr Vec3f Min(Vec3f a, float b)
{
    return Vec3f(Min<float>(a.x, b), Min<float>(a.y, b), Min<float>(a.z, b));
}

YBI_CPU_GPU constexpr Vec3f Min(float a, Vec3f b)
{
    return Min(b, a);
}

YBI_CPU_GPU constexpr Vec3f Max(Vec3f a, Vec3f b)
{
    return Vec3f(Max<float>(a.x, b.x), Max<float>(a.y, b.y), Max<float>(a.z, b.z));
}

YBI_CPU_GPU constexpr Vec3f Max(Vec3f a, float b)
{
    return Vec3f(Max<float>(a.x, b), Max<float>(a.y, b), Max<float>(a.z, b));
}

YBI_CPU_GPU constexpr Vec3f Max(float a, Vec3f b)
{
    return Max(b, a);
}

YBI_CPU_GPU constexpr Vec3f Abs(Vec3f v)
{
    return Vec3f(Abs(v.x), Abs(v.y), Abs(v.z));
}

YBI_CPU_GPU constexpr Vec3f Clamp(const Vec3f &v, float low, float high)
{
    return Max(Min(v, high), low);
}

YBI_CPU_GPU constexpr Vec3f Pow(const Vec3f &v, float power)
{
    return Vec3f(powf(v.x, power), powf(v.y, power), powf(v.z, power));
}

} // namespace ybi
