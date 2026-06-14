#pragma once

#include "core/vec3f.h"
#include <cmath>
#include <cstddef>

namespace ybi
{

class Vec4f
{
public:
    float x;
    float y;
    float z;
    float w;

    constexpr Vec4f() : x(0.0f), y(0.0f), z(0.0f), w(0.0f) {}
    constexpr explicit Vec4f(float value) : x(value), y(value), z(value), w(value) {}
    constexpr Vec4f(float x, float y, float z, float w) : x(x), y(y), z(z), w(w) {}
    constexpr Vec4f(Vec3f xyz, float w) : x(xyz.x), y(xyz.y), z(xyz.z), w(w) {}

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
        YBI_ASSERT(index < 4);
        return Data()[index];
    }

    const float &operator[](size_t index) const
    {
        YBI_ASSERT(index < 4);
        return Data()[index];
    }

    constexpr Vec3f Xyz() const
    {
        return Vec3f(x, y, z);
    }

    constexpr Vec4f operator+() const
    {
        return *this;
    }

    constexpr Vec4f operator-() const
    {
        return Vec4f(-x, -y, -z, -w);
    }

    Vec4f &operator+=(Vec4f rhs)
    {
        x += rhs.x;
        y += rhs.y;
        z += rhs.z;
        w += rhs.w;
        return *this;
    }

    Vec4f &operator-=(Vec4f rhs)
    {
        x -= rhs.x;
        y -= rhs.y;
        z -= rhs.z;
        w -= rhs.w;
        return *this;
    }

    Vec4f &operator*=(float scalar)
    {
        x *= scalar;
        y *= scalar;
        z *= scalar;
        w *= scalar;
        return *this;
    }

    Vec4f &operator/=(float scalar)
    {
        YBI_ASSERT(scalar != 0.0f);
        const float invScalar = 1.0f / scalar;
        return *this *= invScalar;
    }
};

constexpr Vec4f operator+(Vec4f lhs, Vec4f rhs)
{
    return Vec4f(lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z, lhs.w + rhs.w);
}

constexpr Vec4f operator-(Vec4f lhs, Vec4f rhs)
{
    return Vec4f(lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z, lhs.w - rhs.w);
}

constexpr Vec4f operator*(Vec4f lhs, Vec4f rhs)
{
    return Vec4f(lhs.x * rhs.x, lhs.y * rhs.y, lhs.z * rhs.z, lhs.w * rhs.w);
}

constexpr Vec4f operator*(Vec4f vector, float scalar)
{
    return Vec4f(vector.x * scalar, vector.y * scalar, vector.z * scalar, vector.w * scalar);
}

constexpr Vec4f operator*(float scalar, Vec4f vector)
{
    return vector * scalar;
}

inline Vec4f operator/(Vec4f vector, float scalar)
{
    YBI_ASSERT(scalar != 0.0f);
    const float invScalar = 1.0f / scalar;
    return vector * invScalar;
}

constexpr bool operator==(Vec4f lhs, Vec4f rhs)
{
    return lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z && lhs.w == rhs.w;
}

constexpr bool operator!=(Vec4f lhs, Vec4f rhs)
{
    return !(lhs == rhs);
}

constexpr float Dot(Vec4f lhs, Vec4f rhs)
{
    return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z + lhs.w * rhs.w;
}

constexpr float LengthSquared(Vec4f vector)
{
    return Dot(vector, vector);
}

inline float Length(Vec4f vector)
{
    return std::sqrt(LengthSquared(vector));
}

inline Vec4f Normalize(Vec4f vector)
{
    const float length = Length(vector);
    YBI_ASSERT(length != 0.0f);
    return vector / length;
}

constexpr Vec4f Lerp(Vec4f a, Vec4f b, float t)
{
    return a + (b - a) * t;
}

} // namespace ybi
