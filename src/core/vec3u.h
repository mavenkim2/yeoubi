#pragma once

#include "core/assert.h"
#include "core/math.h"
#include "device/cuda_compat.h"
#include <cstddef>
#include <cstdint>

namespace ybi
{

class Vec3u
{
public:
    uint32_t x;
    uint32_t y;
    uint32_t z;

    Vec3u() = default;
    YBI_CPU_GPU constexpr explicit Vec3u(uint32_t value) : x(value), y(value), z(value) {}
    YBI_CPU_GPU constexpr Vec3u(uint32_t x, uint32_t y, uint32_t z) : x(x), y(y), z(z) {}

    YBI_CPU_GPU uint32_t *Data()
    {
        return &x;
    }

    YBI_CPU_GPU const uint32_t *Data() const
    {
        return &x;
    }

    YBI_CPU_GPU uint32_t &operator[](size_t index)
    {
        YBI_ASSERT(index < 3);
        return Data()[index];
    }

    YBI_CPU_GPU const uint32_t &operator[](size_t index) const
    {
        YBI_ASSERT(index < 3);
        return Data()[index];
    }

    YBI_CPU_GPU constexpr Vec3u operator+() const
    {
        return *this;
    }

    YBI_CPU_GPU constexpr Vec3u operator~() const
    {
        return Vec3u(~x, ~y, ~z);
    }

    YBI_CPU_GPU Vec3u &operator+=(Vec3u rhs)
    {
        x += rhs.x;
        y += rhs.y;
        z += rhs.z;
        return *this;
    }

    YBI_CPU_GPU Vec3u &operator-=(Vec3u rhs)
    {
        x -= rhs.x;
        y -= rhs.y;
        z -= rhs.z;
        return *this;
    }

    YBI_CPU_GPU Vec3u &operator*=(uint32_t scalar)
    {
        x *= scalar;
        y *= scalar;
        z *= scalar;
        return *this;
    }

    YBI_CPU_GPU Vec3u &operator/=(uint32_t scalar)
    {
        YBI_ASSERT(scalar != 0);
        x /= scalar;
        y /= scalar;
        z /= scalar;
        return *this;
    }

    YBI_CPU_GPU Vec3u &operator/=(Vec3u rhs)
    {
        YBI_ASSERT(rhs.x != 0);
        YBI_ASSERT(rhs.y != 0);
        YBI_ASSERT(rhs.z != 0);
        x /= rhs.x;
        y /= rhs.y;
        z /= rhs.z;
        return *this;
    }

    YBI_CPU_GPU Vec3u &operator&=(Vec3u rhs)
    {
        x &= rhs.x;
        y &= rhs.y;
        z &= rhs.z;
        return *this;
    }

    YBI_CPU_GPU Vec3u &operator|=(Vec3u rhs)
    {
        x |= rhs.x;
        y |= rhs.y;
        z |= rhs.z;
        return *this;
    }

    YBI_CPU_GPU Vec3u &operator^=(Vec3u rhs)
    {
        x ^= rhs.x;
        y ^= rhs.y;
        z ^= rhs.z;
        return *this;
    }

    YBI_CPU_GPU Vec3u &operator<<=(uint32_t shift)
    {
        x <<= shift;
        y <<= shift;
        z <<= shift;
        return *this;
    }

    YBI_CPU_GPU Vec3u &operator>>=(uint32_t shift)
    {
        x >>= shift;
        y >>= shift;
        z >>= shift;
        return *this;
    }
};

YBI_CPU_GPU constexpr Vec3u operator+(Vec3u lhs, Vec3u rhs)
{
    return Vec3u(lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z);
}

YBI_CPU_GPU constexpr Vec3u operator-(Vec3u lhs, Vec3u rhs)
{
    return Vec3u(lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z);
}

YBI_CPU_GPU constexpr Vec3u operator-(Vec3u lhs, uint32_t rhs)
{
    return Vec3u(lhs.x - rhs, lhs.y - rhs, lhs.z - rhs);
}

YBI_CPU_GPU constexpr Vec3u operator*(Vec3u lhs, Vec3u rhs)
{
    return Vec3u(lhs.x * rhs.x, lhs.y * rhs.y, lhs.z * rhs.z);
}

YBI_CPU_GPU constexpr Vec3u operator*(Vec3u vector, uint32_t scalar)
{
    return Vec3u(vector.x * scalar, vector.y * scalar, vector.z * scalar);
}

YBI_CPU_GPU constexpr Vec3u operator*(uint32_t scalar, Vec3u vector)
{
    return vector * scalar;
}

YBI_CPU_GPU Vec3u operator/(Vec3u vector, uint32_t scalar)
{
    vector /= scalar;
    return vector;
}

YBI_CPU_GPU Vec3u operator/(Vec3u lhs, Vec3u rhs)
{
    lhs /= rhs;
    return lhs;
}

YBI_CPU_GPU constexpr Vec3u operator&(Vec3u lhs, Vec3u rhs)
{
    return Vec3u(lhs.x & rhs.x, lhs.y & rhs.y, lhs.z & rhs.z);
}

YBI_CPU_GPU constexpr Vec3u operator|(Vec3u lhs, Vec3u rhs)
{
    return Vec3u(lhs.x | rhs.x, lhs.y | rhs.y, lhs.z | rhs.z);
}

YBI_CPU_GPU constexpr Vec3u operator^(Vec3u lhs, Vec3u rhs)
{
    return Vec3u(lhs.x ^ rhs.x, lhs.y ^ rhs.y, lhs.z ^ rhs.z);
}

YBI_CPU_GPU constexpr Vec3u operator<<(Vec3u lhs, uint32_t shift)
{
    return Vec3u(lhs.x << shift, lhs.y << shift, lhs.z << shift);
}

YBI_CPU_GPU constexpr Vec3u operator<<(Vec3u lhs, Vec3u rhs)
{
    return Vec3u(lhs.x << rhs.x, lhs.y << rhs.y, lhs.z << rhs.z);
}

YBI_CPU_GPU constexpr Vec3u operator>>(Vec3u lhs, uint32_t shift)
{
    return Vec3u(lhs.x >> shift, lhs.y >> shift, lhs.z >> shift);
}

YBI_CPU_GPU constexpr Vec3u operator>>(Vec3u lhs, Vec3u rhs)
{
    return Vec3u(lhs.x >> rhs.x, lhs.y >> rhs.y, lhs.z >> rhs.z);
}

YBI_CPU_GPU constexpr bool operator==(Vec3u lhs, Vec3u rhs)
{
    return lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z;
}

YBI_CPU_GPU constexpr bool operator!=(Vec3u lhs, Vec3u rhs)
{
    return !(lhs == rhs);
}

YBI_CPU_GPU constexpr Vec3u Min(Vec3u a, Vec3u b)
{
    return Vec3u(ybi::Min<uint32_t>(a.x, b.x),
                 ybi::Min<uint32_t>(a.y, b.y),
                 ybi::Min<uint32_t>(a.z, b.z));
}

YBI_CPU_GPU constexpr Vec3u Max(Vec3u a, Vec3u b)
{
    return Vec3u(ybi::Max<uint32_t>(a.x, b.x),
                 ybi::Max<uint32_t>(a.y, b.y),
                 ybi::Max<uint32_t>(a.z, b.z));
}

} // namespace ybi
