#pragma once

#include "device/util.h"

#if defined(__CUDACC__)
#include <cuda_fp16.h>
#endif

#include <cstdint>

namespace ybi
{
namespace util
{

YBI_DEVICE uint32_t FloatToBits(float value)
{
    union
    {
        float f;
        uint32_t u;
    } bits = {value};
    return bits.u;
}

YBI_DEVICE float BitsToFloat(uint32_t bits)
{
    union
    {
        uint32_t u;
        float f;
    } value = {bits};
    return value.f;
}

YBI_DEVICE uint16_t FloatToHalfBits(float value)
{
#if defined(__CUDA_ARCH__)
    return __half_as_ushort(__float2half(value));
#else
    const uint32_t bits = FloatToBits(value);
    const uint32_t sign = (bits >> 16u) & 0x8000u;
    const uint32_t exponent = (bits >> 23u) & 0xffu;
    uint32_t mantissa = bits & 0x007fffffu;

    if (exponent == 0xffu)
    {
        if (mantissa == 0u)
        {
            return static_cast<uint16_t>(sign | 0x7c00u);
        }
        mantissa >>= 13u;
        return static_cast<uint16_t>(sign | 0x7c00u | mantissa | (mantissa == 0u ? 1u : 0u));
    }

    int32_t halfExponent = static_cast<int32_t>(exponent) - 127 + 15;
    if (halfExponent <= 0)
    {
        if (halfExponent < -10)
        {
            return static_cast<uint16_t>(sign);
        }

        mantissa |= 0x00800000u;
        const uint32_t shift = static_cast<uint32_t>(14 - halfExponent);
        uint32_t halfMantissa = mantissa >> shift;
        const uint32_t roundBit = 1u << (shift - 1u);
        if ((mantissa & roundBit) != 0u &&
            (((mantissa & (roundBit - 1u)) != 0u) || (halfMantissa & 1u) != 0u))
        {
            ++halfMantissa;
        }
        return static_cast<uint16_t>(sign | halfMantissa);
    }

    if (halfExponent >= 31)
    {
        return static_cast<uint16_t>(sign | 0x7c00u);
    }

    uint32_t halfMantissa = mantissa >> 13u;
    if ((mantissa & 0x00001000u) != 0u)
    {
        ++halfMantissa;
        if (halfMantissa == 0x0400u)
        {
            halfMantissa = 0u;
            ++halfExponent;
            if (halfExponent >= 31)
            {
                return static_cast<uint16_t>(sign | 0x7c00u);
            }
        }
    }

    return static_cast<uint16_t>(sign | (static_cast<uint32_t>(halfExponent) << 10u) |
                                 (halfMantissa & 0x03ffu));
#endif
}

YBI_DEVICE float HalfBitsToFloat(uint16_t bits)
{
#if defined(__CUDA_ARCH__)
    return __half2float(__ushort_as_half(bits));
#else
    const uint32_t sign = (static_cast<uint32_t>(bits & 0x8000u)) << 16u;
    const uint32_t exponent = (bits >> 10u) & 0x1fu;
    uint32_t mantissa = bits & 0x03ffu;

    uint32_t floatBits = sign;
    if (exponent == 0u)
    {
        if (mantissa != 0u)
        {
            int32_t normalizedExponent = -14;
            while ((mantissa & 0x0400u) == 0u)
            {
                mantissa <<= 1u;
                --normalizedExponent;
            }
            mantissa &= 0x03ffu;
            floatBits |= static_cast<uint32_t>(normalizedExponent + 127) << 23u;
            floatBits |= mantissa << 13u;
        }
    }
    else if (exponent == 0x1fu)
    {
        floatBits |= 0x7f800000u | (mantissa << 13u);
    }
    else
    {
        floatBits |= (exponent + 112u) << 23u;
        floatBits |= mantissa << 13u;
    }

    return BitsToFloat(floatBits);
#endif
}

} // namespace util
} // namespace ybi
