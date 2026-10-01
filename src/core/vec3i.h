#pragma once

#include "core/assert.h"
#include "core/math.h"
#include "device/cuda_compat.h"

namespace ybi
{

class Vec3i
{
public:
    int32_t x;
    int32_t y;
    int32_t z;

    Vec3i() = default;
    YBI_CPU_GPU constexpr explicit Vec3i(int32_t value) : x(value), y(value), z(value) {}
    YBI_CPU_GPU constexpr Vec3i(int32_t x, int32_t y, int32_t z) : x(x), y(y), z(z) {}

    YBI_CPU_GPU int32_t *Data()
    {
        return &x;
    }

    YBI_CPU_GPU const int32_t *Data() const
    {
        return &x;
    }

    YBI_CPU_GPU int32_t &operator[](size_t index)
    {
        YBI_ASSERT(index < 3);
        return Data()[index];
    }

    YBI_CPU_GPU const int32_t &operator[](size_t index) const
    {
        YBI_ASSERT(index < 3);
        return Data()[index];
    }

    YBI_CPU_GPU constexpr Vec3i operator+() const
    {
        return *this;
    }

    YBI_CPU_GPU constexpr Vec3i operator-() const
    {
        return Vec3i(-x, -y, -z);
    }

    YBI_CPU_GPU constexpr Vec3i operator~() const
    {
        return Vec3i(~x, ~y, ~z);
    }

    YBI_CPU_GPU Vec3i &operator+=(Vec3i rhs)
    {
        x += rhs.x;
        y += rhs.y;
        z += rhs.z;
        return *this;
    }

    YBI_CPU_GPU Vec3i &operator-=(Vec3i rhs)
    {
        x -= rhs.x;
        y -= rhs.y;
        z -= rhs.z;
        return *this;
    }

    YBI_CPU_GPU Vec3i &operator*=(int32_t scalar)
    {
        x *= scalar;
        y *= scalar;
        z *= scalar;
        return *this;
    }

    YBI_CPU_GPU Vec3i &operator/=(int32_t scalar)
    {
        YBI_ASSERT(scalar != 0);
        x /= scalar;
        y /= scalar;
        z /= scalar;
        return *this;
    }

    YBI_CPU_GPU Vec3i &operator/=(Vec3i rhs)
    {
        YBI_ASSERT(rhs.x != 0);
        YBI_ASSERT(rhs.y != 0);
        YBI_ASSERT(rhs.z != 0);
        x /= rhs.x;
        y /= rhs.y;
        z /= rhs.z;
        return *this;
    }

    YBI_CPU_GPU Vec3i &operator&=(Vec3i rhs)
    {
        x &= rhs.x;
        y &= rhs.y;
        z &= rhs.z;
        return *this;
    }

    YBI_CPU_GPU Vec3i &operator|=(Vec3i rhs)
    {
        x |= rhs.x;
        y |= rhs.y;
        z |= rhs.z;
        return *this;
    }

    YBI_CPU_GPU Vec3i &operator^=(Vec3i rhs)
    {
        x ^= rhs.x;
        y ^= rhs.y;
        z ^= rhs.z;
        return *this;
    }

    YBI_CPU_GPU Vec3i &operator<<=(uint32_t shift)
    {
        x <<= shift;
        y <<= shift;
        z <<= shift;
        return *this;
    }

    YBI_CPU_GPU Vec3i &operator>>=(uint32_t shift)
    {
        x >>= shift;
        y >>= shift;
        z >>= shift;
        return *this;
    }
};

YBI_CPU_GPU constexpr Vec3i operator+(Vec3i lhs, Vec3i rhs)
{
    return Vec3i(lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z);
}

YBI_CPU_GPU constexpr Vec3i operator-(Vec3i lhs, Vec3i rhs)
{
    return Vec3i(lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z);
}

YBI_CPU_GPU constexpr Vec3i operator*(Vec3i lhs, Vec3i rhs)
{
    return Vec3i(lhs.x * rhs.x, lhs.y * rhs.y, lhs.z * rhs.z);
}

YBI_CPU_GPU constexpr Vec3i operator*(Vec3i vector, int32_t scalar)
{
    return Vec3i(vector.x * scalar, vector.y * scalar, vector.z * scalar);
}

YBI_CPU_GPU constexpr Vec3i operator*(int32_t scalar, Vec3i vector)
{
    return vector * scalar;
}

YBI_CPU_GPU Vec3i operator/(Vec3i vector, int32_t scalar)
{
    vector /= scalar;
    return vector;
}

YBI_CPU_GPU Vec3i operator/(Vec3i lhs, Vec3i rhs)
{
    lhs /= rhs;
    return lhs;
}

YBI_CPU_GPU constexpr Vec3i operator&(Vec3i lhs, Vec3i rhs)
{
    return Vec3i(lhs.x & rhs.x, lhs.y & rhs.y, lhs.z & rhs.z);
}

YBI_CPU_GPU constexpr Vec3i operator|(Vec3i lhs, Vec3i rhs)
{
    return Vec3i(lhs.x | rhs.x, lhs.y | rhs.y, lhs.z | rhs.z);
}

YBI_CPU_GPU constexpr Vec3i operator^(Vec3i lhs, Vec3i rhs)
{
    return Vec3i(lhs.x ^ rhs.x, lhs.y ^ rhs.y, lhs.z ^ rhs.z);
}

YBI_CPU_GPU constexpr Vec3i operator<<(Vec3i lhs, uint32_t shift)
{
    return Vec3i(lhs.x << shift, lhs.y << shift, lhs.z << shift);
}

YBI_CPU_GPU constexpr Vec3i operator<<(Vec3i lhs, Vec3i rhs)
{
    return Vec3i(lhs.x << rhs.x, lhs.y << rhs.y, lhs.z << rhs.z);
}

YBI_CPU_GPU constexpr Vec3i operator>>(Vec3i lhs, uint32_t shift)
{
    return Vec3i(lhs.x >> shift, lhs.y >> shift, lhs.z >> shift);
}

YBI_CPU_GPU constexpr Vec3i operator>>(Vec3i lhs, Vec3i rhs)
{
    return Vec3i(lhs.x >> rhs.x, lhs.y >> rhs.y, lhs.z >> rhs.z);
}

YBI_CPU_GPU constexpr bool operator==(Vec3i lhs, Vec3i rhs)
{
    return lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z;
}

YBI_CPU_GPU constexpr bool operator!=(Vec3i lhs, Vec3i rhs)
{
    return !(lhs == rhs);
}

YBI_CPU_GPU constexpr Vec3i Min(Vec3i a, Vec3i b)
{
    return Vec3i(ybi::Min<int32_t>(a.x, b.x),
                 ybi::Min<int32_t>(a.y, b.y),
                 ybi::Min<int32_t>(a.z, b.z));
}

YBI_CPU_GPU constexpr Vec3i Max(Vec3i a, Vec3i b)
{
    return Vec3i(ybi::Max<int32_t>(a.x, b.x),
                 ybi::Max<int32_t>(a.y, b.y),
                 ybi::Max<int32_t>(a.z, b.z));
}

} // namespace ybi
