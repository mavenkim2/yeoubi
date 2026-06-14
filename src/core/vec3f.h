#pragma once

#include "core/assert.h"
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

    constexpr Vec3f() : x(0.0f), y(0.0f), z(0.0f) {}
    constexpr explicit Vec3f(float value) : x(value), y(value), z(value) {}
    constexpr Vec3f(float x, float y, float z) : x(x), y(y), z(z) {}

    float *Data()
    {
        return &x;
    }

    const float *Data() const
    {
        return &x;
    }

    float &operator[](size_t index)
    {
        YBI_ASSERT(index < 3);
        return Data()[index];
    }

    const float &operator[](size_t index) const
    {
        YBI_ASSERT(index < 3);
        return Data()[index];
    }

    constexpr Vec3f operator+() const
    {
        return *this;
    }

    constexpr Vec3f operator-() const
    {
        return Vec3f(-x, -y, -z);
    }

    Vec3f &operator+=(Vec3f rhs)
    {
        x += rhs.x;
        y += rhs.y;
        z += rhs.z;
        return *this;
    }

    Vec3f &operator-=(Vec3f rhs)
    {
        x -= rhs.x;
        y -= rhs.y;
        z -= rhs.z;
        return *this;
    }

    Vec3f &operator*=(float scalar)
    {
        x *= scalar;
        y *= scalar;
        z *= scalar;
        return *this;
    }

    Vec3f &operator/=(float scalar)
    {
        YBI_ASSERT(scalar != 0.0f);
        const float invScalar = 1.0f / scalar;
        return *this *= invScalar;
    }
};

constexpr Vec3f operator+(Vec3f lhs, Vec3f rhs)
{
    return Vec3f(lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z);
}

constexpr Vec3f operator-(Vec3f lhs, Vec3f rhs)
{
    return Vec3f(lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z);
}

constexpr Vec3f operator*(Vec3f lhs, Vec3f rhs)
{
    return Vec3f(lhs.x * rhs.x, lhs.y * rhs.y, lhs.z * rhs.z);
}

constexpr Vec3f operator*(Vec3f vector, float scalar)
{
    return Vec3f(vector.x * scalar, vector.y * scalar, vector.z * scalar);
}

constexpr Vec3f operator*(float scalar, Vec3f vector)
{
    return vector * scalar;
}

inline Vec3f operator/(Vec3f vector, float scalar)
{
    YBI_ASSERT(scalar != 0.0f);
    const float invScalar = 1.0f / scalar;
    return vector * invScalar;
}

constexpr bool operator==(Vec3f lhs, Vec3f rhs)
{
    return lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z;
}

constexpr bool operator!=(Vec3f lhs, Vec3f rhs)
{
    return !(lhs == rhs);
}

constexpr float Dot(Vec3f lhs, Vec3f rhs)
{
    return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z;
}

constexpr Vec3f Cross(Vec3f lhs, Vec3f rhs)
{
    return Vec3f(lhs.y * rhs.z - lhs.z * rhs.y,
                 lhs.z * rhs.x - lhs.x * rhs.z,
                 lhs.x * rhs.y - lhs.y * rhs.x);
}

constexpr float LengthSquared(Vec3f vector)
{
    return Dot(vector, vector);
}

inline float Length(Vec3f vector)
{
    return std::sqrt(LengthSquared(vector));
}

inline Vec3f Normalize(Vec3f vector)
{
    const float length = Length(vector);
    YBI_ASSERT(length != 0.0f);
    return vector / length;
}

constexpr Vec3f Lerp(Vec3f a, Vec3f b, float t)
{
    return a + (b - a) * t;
}

} // namespace ybi
