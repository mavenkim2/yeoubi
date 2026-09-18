#pragma once

#include "core/vec2f.h"
#include "core/vec3f.h"
#include "device/cuda_compat.h"
#include <cmath>
#include <cstddef>

namespace ybi
{

class Vec4f
{
public:
    union
    {
        struct
        {
            float x;
            float y;
            float z;
            float w;
        };
        Vec2f xy;
        Vec3f xyz;
    };

    YBI_CPU_GPU constexpr Vec4f() : x(0.0f), y(0.0f), z(0.0f), w(0.0f) {}
    YBI_CPU_GPU constexpr explicit Vec4f(float value) : x(value), y(value), z(value), w(value) {}
    YBI_CPU_GPU constexpr Vec4f(float x, float y, float z, float w) : x(x), y(y), z(z), w(w) {}
    YBI_CPU_GPU constexpr Vec4f(Vec3f xyz, float w) : x(xyz.x), y(xyz.y), z(xyz.z), w(w) {}

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
        YBI_ASSERT(index < 4);
        return Data()[index];
    }

    YBI_CPU_GPU const float &operator[](size_t index) const
    {
        YBI_ASSERT(index < 4);
        return Data()[index];
    }

    YBI_CPU_GPU constexpr Vec2f Xy() const
    {
        return Vec2f(x, y);
    }

    YBI_CPU_GPU constexpr Vec3f Xyz() const
    {
        return Vec3f(x, y, z);
    }

    YBI_CPU_GPU constexpr Vec2f Zw() const
    {
        return Vec2f(z, w);
    }

    YBI_CPU_GPU constexpr Vec4f operator+() const
    {
        return *this;
    }

    YBI_CPU_GPU constexpr Vec4f operator-() const
    {
        return Vec4f(-x, -y, -z, -w);
    }

    YBI_CPU_GPU Vec4f &operator+=(Vec4f rhs)
    {
        x += rhs.x;
        y += rhs.y;
        z += rhs.z;
        w += rhs.w;
        return *this;
    }

    YBI_CPU_GPU Vec4f &operator-=(Vec4f rhs)
    {
        x -= rhs.x;
        y -= rhs.y;
        z -= rhs.z;
        w -= rhs.w;
        return *this;
    }

    YBI_CPU_GPU Vec4f &operator*=(float scalar)
    {
        x *= scalar;
        y *= scalar;
        z *= scalar;
        w *= scalar;
        return *this;
    }

    YBI_CPU_GPU Vec4f &operator/=(float scalar)
    {
        YBI_ASSERT(scalar != 0.0f);
        const float invScalar = 1.0f / scalar;
        return *this *= invScalar;
    }
};

YBI_CPU_GPU constexpr Vec4f operator+(Vec4f lhs, Vec4f rhs)
{
    return Vec4f(lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z, lhs.w + rhs.w);
}

YBI_CPU_GPU constexpr Vec4f operator-(Vec4f lhs, Vec4f rhs)
{
    return Vec4f(lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z, lhs.w - rhs.w);
}

YBI_CPU_GPU constexpr Vec4f operator*(Vec4f lhs, Vec4f rhs)
{
    return Vec4f(lhs.x * rhs.x, lhs.y * rhs.y, lhs.z * rhs.z, lhs.w * rhs.w);
}

YBI_CPU_GPU constexpr Vec4f operator*(Vec4f vector, float scalar)
{
    return Vec4f(vector.x * scalar, vector.y * scalar, vector.z * scalar, vector.w * scalar);
}

YBI_CPU_GPU constexpr Vec4f operator*(float scalar, Vec4f vector)
{
    return vector * scalar;
}

YBI_CPU_GPU Vec4f operator/(Vec4f vector, float scalar)
{
    YBI_ASSERT(scalar != 0.0f);
    const float invScalar = 1.0f / scalar;
    return vector * invScalar;
}

YBI_CPU_GPU constexpr bool operator==(Vec4f lhs, Vec4f rhs)
{
    return lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z && lhs.w == rhs.w;
}

YBI_CPU_GPU constexpr bool operator!=(Vec4f lhs, Vec4f rhs)
{
    return !(lhs == rhs);
}

YBI_CPU_GPU constexpr float Dot(Vec4f lhs, Vec4f rhs)
{
    return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z + lhs.w * rhs.w;
}

YBI_CPU_GPU constexpr float LengthSquared(Vec4f vector)
{
    return Dot(vector, vector);
}

YBI_CPU_GPU float Length(Vec4f vector)
{
    return std::sqrt(LengthSquared(vector));
}

YBI_CPU_GPU Vec4f Normalize(Vec4f vector)
{
    const float length = Length(vector);
    YBI_ASSERT(length != 0.0f);
    return vector / length;
}

YBI_CPU_GPU constexpr Vec4f Lerp(Vec4f a, Vec4f b, float t)
{
    return a + (b - a) * t;
}

} // namespace ybi
