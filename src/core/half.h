#pragma once

#include <cstdint>
#include <cstring>
#include <limits>

namespace ybi
{

// IEEE 754 binary16 value. Arithmetic is performed as binary32 on the host and
// rounded back to binary16 using round-to-nearest, ties-to-even.
class Half
{
public:
    Half() = default;
    Half(float value) : bits(FloatToBits(value)) {}

    static Half FromBits(uint16_t value)
    {
        Half result;
        result.bits = value;
        return result;
    }

    uint16_t Bits() const
    {
        return bits;
    }

    operator float() const
    {
        return BitsToFloat(bits);
    }

    Half &operator+=(Half rhs)
    {
        *this = Half(static_cast<float>(*this) + static_cast<float>(rhs));
        return *this;
    }

    Half &operator-=(Half rhs)
    {
        *this = Half(static_cast<float>(*this) - static_cast<float>(rhs));
        return *this;
    }

    Half &operator*=(Half rhs)
    {
        *this = Half(static_cast<float>(*this) * static_cast<float>(rhs));
        return *this;
    }

    Half &operator/=(Half rhs)
    {
        *this = Half(static_cast<float>(*this) / static_cast<float>(rhs));
        return *this;
    }

    friend Half operator+(Half lhs, Half rhs)
    {
        lhs += rhs;
        return lhs;
    }

    friend Half operator-(Half lhs, Half rhs)
    {
        lhs -= rhs;
        return lhs;
    }

    friend Half operator*(Half lhs, Half rhs)
    {
        lhs *= rhs;
        return lhs;
    }

    friend Half operator/(Half lhs, Half rhs)
    {
        lhs /= rhs;
        return lhs;
    }

    Half operator-() const
    {
        return FromBits(static_cast<uint16_t>(bits ^ 0x8000u));
    }

private:
    static uint32_t FloatBits(float value)
    {
        uint32_t result;
        static_assert(sizeof(result) == sizeof(value));
        std::memcpy(&result, &value, sizeof(result));
        return result;
    }

    static float FloatFromBits(uint32_t value)
    {
        float result;
        static_assert(sizeof(result) == sizeof(value));
        std::memcpy(&result, &value, sizeof(result));
        return result;
    }

    static uint16_t FloatToBits(float value)
    {
        const uint32_t floatBits = FloatBits(value);
        const uint16_t sign = static_cast<uint16_t>((floatBits >> 16u) & 0x8000u);
        const uint32_t exponent = (floatBits >> 23u) & 0xffu;
        const uint32_t mantissa = floatBits & 0x7fffffu;

        if (exponent == 0xffu)
        {
            if (mantissa == 0)
            {
                return static_cast<uint16_t>(sign | 0x7c00u);
            }

            // Preserve as much of the NaN payload as binary16 can hold. Ensure
            // it remains a NaN if all retained payload bits were zero.
            uint16_t payload = static_cast<uint16_t>(mantissa >> 13u);
            if (payload == 0)
            {
                payload = 1;
            }
            return static_cast<uint16_t>(sign | 0x7c00u | payload);
        }

        const int halfExponent = static_cast<int>(exponent) - 127 + 15;
        if (halfExponent >= 31)
        {
            return static_cast<uint16_t>(sign | 0x7c00u);
        }

        if (halfExponent <= 0)
        {
            if (halfExponent < -10)
            {
                return sign;
            }

            const uint32_t significand = mantissa | 0x800000u;
            const uint32_t shift = static_cast<uint32_t>(14 - halfExponent);
            uint32_t halfMantissa = significand >> shift;
            const uint32_t remainderMask = (uint32_t{1} << shift) - 1u;
            const uint32_t remainder = significand & remainderMask;
            const uint32_t halfway = uint32_t{1} << (shift - 1u);

            if (remainder > halfway || (remainder == halfway && (halfMantissa & 1u) != 0))
            {
                ++halfMantissa;
            }

            return static_cast<uint16_t>(sign | halfMantissa);
        }

        uint32_t result = static_cast<uint32_t>(sign) |
                          (static_cast<uint32_t>(halfExponent) << 10u) | (mantissa >> 13u);
        const uint32_t remainder = mantissa & 0x1fffu;
        if (remainder > 0x1000u || (remainder == 0x1000u && (result & 1u) != 0))
        {
            ++result;
        }

        return static_cast<uint16_t>(result);
    }

    static float BitsToFloat(uint16_t value)
    {
        const uint32_t sign = (static_cast<uint32_t>(value) & 0x8000u) << 16u;
        const uint32_t exponent = (static_cast<uint32_t>(value) >> 10u) & 0x1fu;
        uint32_t mantissa = static_cast<uint32_t>(value) & 0x3ffu;
        uint32_t result;

        if (exponent == 0)
        {
            if (mantissa == 0)
            {
                result = sign;
            }
            else
            {
                int normalizedExponent = -14;
                while ((mantissa & 0x400u) == 0)
                {
                    mantissa <<= 1u;
                    --normalizedExponent;
                }
                mantissa &= 0x3ffu;
                result = sign | (static_cast<uint32_t>(normalizedExponent + 127) << 23u) |
                         (mantissa << 13u);
            }
        }
        else if (exponent == 0x1fu)
        {
            result = sign | 0x7f800000u | (mantissa << 13u);
        }
        else
        {
            result = sign | ((exponent + (127u - 15u)) << 23u) | (mantissa << 13u);
        }

        return FloatFromBits(result);
    }

    uint16_t bits = 0;
};

static_assert(sizeof(Half) == sizeof(uint16_t));
static_assert(sizeof(float) == sizeof(uint32_t));
static_assert(std::numeric_limits<float>::is_iec559);

} // namespace ybi
