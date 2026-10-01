#pragma once

#include "core/assert.h"
#include "core/math.h"
#include "device/cuda_compat.h"

namespace ybi
{

class Vec4u
{
public:
    uint32_t x;
    uint32_t y;
    uint32_t z;
    uint32_t w;

    Vec4u() = default;
    YBI_CPU_GPU constexpr explicit Vec4u(uint32_t value) : x(value), y(value), z(value), w(value) {}
    YBI_CPU_GPU constexpr Vec4u(uint32_t x, uint32_t y, uint32_t z, uint32_t w)
        : x(x), y(y), z(z), w(w)
    {
    }

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
        YBI_ASSERT(index < 4);
        return Data()[index];
    }

    YBI_CPU_GPU const uint32_t &operator[](size_t index) const
    {
        YBI_ASSERT(index < 4);
        return Data()[index];
    }

    YBI_CPU_GPU constexpr Vec4u operator+() const
    {
        return *this;
    }

    YBI_CPU_GPU constexpr Vec4u operator~() const
    {
        return Vec4u(~x, ~y, ~z, ~w);
    }

    YBI_CPU_GPU Vec4u &operator+=(Vec4u rhs)
    {
        x += rhs.x;
        y += rhs.y;
        z += rhs.z;
        w += rhs.w;
        return *this;
    }

    YBI_CPU_GPU Vec4u &operator-=(Vec4u rhs)
    {
        x -= rhs.x;
        y -= rhs.y;
        z -= rhs.z;
        w -= rhs.w;
        return *this;
    }

    YBI_CPU_GPU Vec4u &operator*=(uint32_t scalar)
    {
        x *= scalar;
        y *= scalar;
        z *= scalar;
        w *= scalar;
        return *this;
    }

    YBI_CPU_GPU Vec4u &operator/=(uint32_t scalar)
    {
        YBI_ASSERT(scalar != 0);
        x /= scalar;
        y /= scalar;
        z /= scalar;
        w /= scalar;
        return *this;
    }

    YBI_CPU_GPU Vec4u &operator/=(Vec4u rhs)
    {
        YBI_ASSERT(rhs.x != 0);
        YBI_ASSERT(rhs.y != 0);
        YBI_ASSERT(rhs.z != 0);
        YBI_ASSERT(rhs.w != 0);
        x /= rhs.x;
        y /= rhs.y;
        z /= rhs.z;
        w /= rhs.w;
        return *this;
    }

    YBI_CPU_GPU Vec4u &operator&=(Vec4u rhs)
    {
        x &= rhs.x;
        y &= rhs.y;
        z &= rhs.z;
        w &= rhs.w;
        return *this;
    }

    YBI_CPU_GPU Vec4u &operator|=(Vec4u rhs)
    {
        x |= rhs.x;
        y |= rhs.y;
        z |= rhs.z;
        w |= rhs.w;
        return *this;
    }

    YBI_CPU_GPU Vec4u &operator^=(Vec4u rhs)
    {
        x ^= rhs.x;
        y ^= rhs.y;
        z ^= rhs.z;
        w ^= rhs.w;
        return *this;
    }

    YBI_CPU_GPU Vec4u &operator<<=(uint32_t shift)
    {
        x <<= shift;
        y <<= shift;
        z <<= shift;
        w <<= shift;
        return *this;
    }

    YBI_CPU_GPU Vec4u &operator>>=(uint32_t shift)
    {
        x >>= shift;
        y >>= shift;
        z >>= shift;
        w >>= shift;
        return *this;
    }
};

YBI_CPU_GPU constexpr Vec4u operator+(Vec4u lhs, Vec4u rhs)
{
    return Vec4u(lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z, lhs.w + rhs.w);
}

YBI_CPU_GPU constexpr Vec4u operator-(Vec4u lhs, Vec4u rhs)
{
    return Vec4u(lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z, lhs.w - rhs.w);
}

YBI_CPU_GPU constexpr Vec4u operator-(Vec4u lhs, uint32_t rhs)
{
    return Vec4u(lhs.x - rhs, lhs.y - rhs, lhs.z - rhs, lhs.w - rhs);
}

YBI_CPU_GPU constexpr Vec4u operator*(Vec4u lhs, Vec4u rhs)
{
    return Vec4u(lhs.x * rhs.x, lhs.y * rhs.y, lhs.z * rhs.z, lhs.w * rhs.w);
}

YBI_CPU_GPU constexpr Vec4u operator*(Vec4u vector, uint32_t scalar)
{
    return Vec4u(vector.x * scalar, vector.y * scalar, vector.z * scalar, vector.w * scalar);
}

YBI_CPU_GPU constexpr Vec4u operator*(uint32_t scalar, Vec4u vector)
{
    return vector * scalar;
}

YBI_CPU_GPU Vec4u operator/(Vec4u vector, uint32_t scalar)
{
    vector /= scalar;
    return vector;
}

YBI_CPU_GPU Vec4u operator/(Vec4u lhs, Vec4u rhs)
{
    lhs /= rhs;
    return lhs;
}

YBI_CPU_GPU constexpr Vec4u operator&(Vec4u lhs, Vec4u rhs)
{
    return Vec4u(lhs.x & rhs.x, lhs.y & rhs.y, lhs.z & rhs.z, lhs.w & rhs.w);
}

YBI_CPU_GPU constexpr Vec4u operator|(Vec4u lhs, Vec4u rhs)
{
    return Vec4u(lhs.x | rhs.x, lhs.y | rhs.y, lhs.z | rhs.z, lhs.w | rhs.w);
}

YBI_CPU_GPU constexpr Vec4u operator^(Vec4u lhs, Vec4u rhs)
{
    return Vec4u(lhs.x ^ rhs.x, lhs.y ^ rhs.y, lhs.z ^ rhs.z, lhs.w ^ rhs.w);
}

YBI_CPU_GPU constexpr Vec4u operator<<(Vec4u lhs, uint32_t shift)
{
    return Vec4u(lhs.x << shift, lhs.y << shift, lhs.z << shift, lhs.w << shift);
}

YBI_CPU_GPU constexpr Vec4u operator<<(Vec4u lhs, Vec4u rhs)
{
    return Vec4u(lhs.x << rhs.x, lhs.y << rhs.y, lhs.z << rhs.z, lhs.w << rhs.w);
}

YBI_CPU_GPU constexpr Vec4u operator>>(Vec4u lhs, uint32_t shift)
{
    return Vec4u(lhs.x >> shift, lhs.y >> shift, lhs.z >> shift, lhs.w >> shift);
}

YBI_CPU_GPU constexpr Vec4u operator>>(Vec4u lhs, Vec4u rhs)
{
    return Vec4u(lhs.x >> rhs.x, lhs.y >> rhs.y, lhs.z >> rhs.z, lhs.w >> rhs.w);
}

YBI_CPU_GPU constexpr bool operator==(Vec4u lhs, Vec4u rhs)
{
    return lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z && lhs.w == rhs.w;
}

YBI_CPU_GPU constexpr bool operator!=(Vec4u lhs, Vec4u rhs)
{
    return !(lhs == rhs);
}

YBI_CPU_GPU constexpr Vec4u Min(Vec4u a, Vec4u b)
{
    return Vec4u(ybi::Min<uint32_t>(a.x, b.x),
                 ybi::Min<uint32_t>(a.y, b.y),
                 ybi::Min<uint32_t>(a.z, b.z),
                 ybi::Min<uint32_t>(a.w, b.w));
}

YBI_CPU_GPU constexpr Vec4u Max(Vec4u a, Vec4u b)
{
    return Vec4u(ybi::Max<uint32_t>(a.x, b.x),
                 ybi::Max<uint32_t>(a.y, b.y),
                 ybi::Max<uint32_t>(a.z, b.z),
                 ybi::Max<uint32_t>(a.w, b.w));
}

} // namespace ybi
