#pragma once

#include "core/assert.h"
#include "core/math.h"
#include "device/cuda_compat.h"

namespace ybi
{

class Vec2f
{
public:
    float x;
    float y;

    YBI_CPU_GPU constexpr Vec2f() : x(0.0f), y(0.0f) {}
    YBI_CPU_GPU constexpr explicit Vec2f(float value) : x(value), y(value) {}
    YBI_CPU_GPU constexpr Vec2f(float x, float y) : x(x), y(y) {}

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
        YBI_ASSERT(index < 2);
        return Data()[index];
    }

    YBI_CPU_GPU const float &operator[](size_t index) const
    {
        YBI_ASSERT(index < 2);
        return Data()[index];
    }

    YBI_CPU_GPU constexpr Vec2f operator+() const
    {
        return *this;
    }

    YBI_CPU_GPU constexpr Vec2f operator-() const
    {
        return Vec2f(-x, -y);
    }

    YBI_CPU_GPU Vec2f &operator+=(Vec2f rhs)
    {
        x += rhs.x;
        y += rhs.y;
        return *this;
    }

    YBI_CPU_GPU Vec2f &operator-=(Vec2f rhs)
    {
        x -= rhs.x;
        y -= rhs.y;
        return *this;
    }

    YBI_CPU_GPU Vec2f &operator*=(float scalar)
    {
        x *= scalar;
        y *= scalar;
        return *this;
    }

    YBI_CPU_GPU Vec2f &operator/=(float scalar)
    {
        YBI_ASSERT(scalar != 0.0f);
        const float invScalar = 1.0f / scalar;
        return *this *= invScalar;
    }
};

YBI_CPU_GPU constexpr Vec2f operator+(Vec2f lhs, Vec2f rhs)
{
    return Vec2f(lhs.x + rhs.x, lhs.y + rhs.y);
}

YBI_CPU_GPU constexpr Vec2f operator+(Vec2f lhs, float rhs)
{
    return Vec2f(lhs.x + rhs, lhs.y + rhs);
}

YBI_CPU_GPU constexpr Vec2f operator+(float lhs, Vec2f rhs)
{
    return rhs + lhs;
}

YBI_CPU_GPU constexpr Vec2f operator-(Vec2f lhs, Vec2f rhs)
{
    return Vec2f(lhs.x - rhs.x, lhs.y - rhs.y);
}

YBI_CPU_GPU constexpr Vec2f operator-(Vec2f lhs, float rhs)
{
    return Vec2f(lhs.x - rhs, lhs.y - rhs);
}

YBI_CPU_GPU constexpr Vec2f operator*(Vec2f lhs, Vec2f rhs)
{
    return Vec2f(lhs.x * rhs.x, lhs.y * rhs.y);
}

YBI_CPU_GPU constexpr Vec2f operator*(Vec2f vector, float scalar)
{
    return Vec2f(vector.x * scalar, vector.y * scalar);
}

YBI_CPU_GPU constexpr Vec2f operator*(float scalar, Vec2f vector)
{
    return vector * scalar;
}

YBI_CPU_GPU Vec2f operator/(Vec2f lhs, Vec2f rhs)
{
    return Vec2f(lhs.x / rhs.x, lhs.y / rhs.y);
}

YBI_CPU_GPU Vec2f operator/(Vec2f vector, float scalar)
{
    YBI_ASSERT(scalar != 0.0f);
    const float invScalar = 1.0f / scalar;
    return vector * invScalar;
}

YBI_CPU_GPU constexpr bool operator==(Vec2f lhs, Vec2f rhs)
{
    return lhs.x == rhs.x && lhs.y == rhs.y;
}

YBI_CPU_GPU constexpr bool operator!=(Vec2f lhs, Vec2f rhs)
{
    return !(lhs == rhs);
}

YBI_CPU_GPU constexpr float Dot(Vec2f lhs, Vec2f rhs)
{
    return lhs.x * rhs.x + lhs.y * rhs.y;
}

YBI_CPU_GPU constexpr float LengthSquared(Vec2f vector)
{
    return Dot(vector, vector);
}

YBI_CPU_GPU float Length(Vec2f vector)
{
    return sqrt(LengthSquared(vector));
}

YBI_CPU_GPU Vec2f Normalize(Vec2f vector)
{
    const float length = Length(vector);
    YBI_ASSERT(length != 0.0f);
    return vector / length;
}

YBI_CPU_GPU constexpr Vec2f Lerp(Vec2f a, Vec2f b, float t)
{
    return a + (b - a) * t;
}

YBI_CPU_GPU Vec2f Floor(Vec2f a)
{
    return Vec2f(floorf(a.x), floorf(a.y));
}

YBI_CPU_GPU constexpr Vec2f Max(const Vec2f &lhs, float rhs)
{
    return Vec2f(ybi::Max<float>(lhs.x, rhs), ybi::Max<float>(lhs.y, rhs));
}

YBI_CPU_GPU constexpr Vec2f Max(float lhs, const Vec2f &rhs)
{
    return Max(rhs, lhs);
}

YBI_CPU_GPU constexpr Vec2f Min(const Vec2f &lhs, float rhs)
{
    return Vec2f(ybi::Min<float>(lhs.x, rhs), ybi::Min<float>(lhs.y, rhs));
}

YBI_CPU_GPU constexpr Vec2f Min(float lhs, const Vec2f &rhs)
{
    return Min(rhs, lhs);
}

YBI_CPU_GPU constexpr Vec2f Clamp(const Vec2f &v, float low, float high)
{
    return Max(Min(v, high), low);
}

} // namespace ybi
