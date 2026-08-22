#pragma once

#include "core/assert.h"
#include "core/math.h"
#include "device/cuda_compat.h"
#include <cstddef>
#include <cstdint>

namespace ybi
{

class Vec2u
{
public:
    uint32_t x;
    uint32_t y;

    Vec2u() = default;
    YBI_CPU_GPU constexpr explicit Vec2u(uint32_t value) : x(value), y(value) {}
    YBI_CPU_GPU constexpr Vec2u(uint32_t x, uint32_t y) : x(x), y(y) {}

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
        YBI_ASSERT(index < 2);
        return Data()[index];
    }

    YBI_CPU_GPU const uint32_t &operator[](size_t index) const
    {
        YBI_ASSERT(index < 2);
        return Data()[index];
    }

    YBI_CPU_GPU constexpr Vec2u operator+() const
    {
        return *this;
    }

    YBI_CPU_GPU Vec2u &operator+=(Vec2u rhs)
    {
        x += rhs.x;
        y += rhs.y;
        return *this;
    }

    YBI_CPU_GPU Vec2u &operator-=(Vec2u rhs)
    {
        x -= rhs.x;
        y -= rhs.y;
        return *this;
    }

    YBI_CPU_GPU Vec2u &operator*=(uint32_t scalar)
    {
        x *= scalar;
        y *= scalar;
        return *this;
    }

    YBI_CPU_GPU Vec2u &operator/=(uint32_t scalar)
    {
        YBI_ASSERT(scalar != 0);
        x /= scalar;
        y /= scalar;
        return *this;
    }

    YBI_CPU_GPU Vec2u &operator^=(uint32_t scalar)
    {
        x ^= scalar;
        y ^= scalar;
        return *this;
    }
};

YBI_CPU_GPU constexpr Vec2u operator+(Vec2u lhs, Vec2u rhs)
{
    return Vec2u(lhs.x + rhs.x, lhs.y + rhs.y);
}

YBI_CPU_GPU constexpr Vec2u operator+(Vec2u lhs, uint32_t rhs)
{
    return Vec2u(lhs.x + rhs, lhs.y + rhs);
}

YBI_CPU_GPU constexpr Vec2u operator+(uint32_t lhs, Vec2u rhs)
{
    return rhs + lhs;
}

YBI_CPU_GPU constexpr Vec2u operator-(Vec2u lhs, Vec2u rhs)
{
    return Vec2u(lhs.x - rhs.x, lhs.y - rhs.y);
}

YBI_CPU_GPU constexpr Vec2u operator*(Vec2u lhs, Vec2u rhs)
{
    return Vec2u(lhs.x * rhs.x, lhs.y * rhs.y);
}

YBI_CPU_GPU constexpr Vec2u operator*(Vec2u vector, uint32_t scalar)
{
    return Vec2u(vector.x * scalar, vector.y * scalar);
}

YBI_CPU_GPU constexpr Vec2u operator*(uint32_t scalar, Vec2u vector)
{
    return vector * scalar;
}

YBI_CPU_GPU Vec2u operator/(Vec2u vector, uint32_t scalar)
{
    vector /= scalar;
    return vector;
}

YBI_CPU_GPU Vec2u operator^(Vec2u lhs, uint32_t rhs)
{
    return Vec2u(lhs.x ^ rhs, lhs.y ^ rhs);
}

YBI_CPU_GPU constexpr bool operator==(Vec2u lhs, Vec2u rhs)
{
    return lhs.x == rhs.x && lhs.y == rhs.y;
}

YBI_CPU_GPU constexpr bool operator!=(Vec2u lhs, Vec2u rhs)
{
    return !(lhs == rhs);
}

YBI_CPU_GPU constexpr Vec2u Min(Vec2u a, Vec2u b)
{
    return Vec2u(ybi::Min<uint32_t>(a.x, b.x), ybi::Min<uint32_t>(a.y, b.y));
}

YBI_CPU_GPU constexpr Vec2u Max(Vec2u a, Vec2u b)
{
    return Vec2u(ybi::Max<uint32_t>(a.x, b.x), ybi::Max<uint32_t>(a.y, b.y));
}

} // namespace ybi
