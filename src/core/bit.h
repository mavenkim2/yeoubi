#pragma once

#include "core/assert.h"
#include "device/cuda_compat.h"
#include <cstdint>
#include <cstring>

#if defined(_MSC_VER) && !defined(__CUDA_ARCH__)
#include <intrin.h>
#endif

namespace ybi
{

// Returns the zero-based index of the most significant set bit, or -1 when value is zero.
YBI_CPU_GPU int FindMostSignificantBit(uint32_t value)
{
    if (value == 0)
    {
        return -1;
    }

#if defined(__CUDA_ARCH__)
    return 31 - __clz(value);
#elif defined(_MSC_VER)
    unsigned long index;
    _BitScanReverse(&index, value);
    return static_cast<int>(index);
#elif defined(__GNUC__) || defined(__clang__)
    return 31 - __builtin_clz(value);
#else
    int index = 0;
    while (value >>= 1)
    {
        index++;
    }
    return index;
#endif
}

YBI_CPU_GPU uint32_t PopCount(uint32_t value)
{
#if defined(__CUDA_ARCH__)
    return __popc(value);
#elif defined(_MSC_VER)
    return static_cast<uint32_t>(__popcnt(value));
#elif defined(__GNUC__) || defined(__clang__)
    return static_cast<uint32_t>(__builtin_popcount(value));
#else
    uint32_t count = 0;
    while (value)
    {
        value &= value - 1u;
        count++;
    }
    return count;
#endif
}

YBI_CPU_GPU bool IsPowerOfTwo(size_t value)
{
    return value != 0 && (value & (value - 1)) == 0;
}

YBI_CPU_GPU size_t AlignUp(size_t val, size_t alignment)
{
    YBI_ASSERT(IsPowerOfTwo(alignment));
    return (val + alignment - 1) & ~(alignment - 1);
}

YBI_CPU_GPU size_t NextPowerOfTwo(size_t v)
{
    v--;
    v |= v >> 1;
    v |= v >> 2;
    v |= v >> 4;
    v |= v >> 8;
    v |= v >> 16;
    v |= v >> 32;
    v++;
    return v;
}

// Returns the low 32 bits of the concatenated high:low value shifted right.
YBI_CPU_GPU uint32_t FunnelShiftRight(uint32_t low, uint32_t high, uint32_t shift)
{
    YBI_ASSERT(shift < 32);

#if defined(__CUDA_ARCH__)
    return __funnelshift_r(low, high, shift);
#else
    return shift == 0 ? low : (low >> shift) | (high << (32 - shift));
#endif
}

YBI_CPU_GPU uint32_t ExtractBits(uint32_t data, uint32_t size, uint32_t offset)
{
    size &= 31;
    offset &= 31;

#if defined(__CUDA_ARCH__)
    uint32_t result;
    asm("bfe.u32 %0, %1, %2, %3;" : "=r"(result) : "r"(data), "r"(offset), "r"(size));
    return result;
#else
    return (data >> offset) & ((1u << size) - 1u);
#endif
}

YBI_CPU_GPU uint32_t PackBits(uint32_t val, uint32_t data, uint32_t &offset, uint32_t size)
{
    uint32_t o = offset & 31u;
    data = data & ((1u << size) - 1u);
    val |= data << o;
    offset += size;
    return val;
}

YBI_CPU_GPU void WriteBits(uint32_t *data, uint32_t &position, uint32_t value, uint32_t numBits)
{
    if (numBits == 0)
        return;
    YBI_ASSERT(numBits <= 32);
    uint32_t dwordIndex = position >> 5;
    uint32_t bitIndex = position & 31;

    YBI_ASSERT(numBits == 32 || ((value & ((1u << numBits) - 1)) == value));

    data[dwordIndex] |= value << bitIndex;
    if (bitIndex + numBits > 32)
    {
        data[dwordIndex + 1] |= value >> (32 - bitIndex);
    }

    position += numBits;
}

YBI_CPU_GPU void WriteBits(uint8_t *data, uint32_t &position, uint32_t value, uint32_t numBits)
{
    if (numBits == 0)
        return;
    YBI_ASSERT(numBits <= 32);
    YBI_ASSERT(numBits == 32 || ((value & ((1u << numBits) - 1)) == value));

    while (numBits)
    {
        const uint32_t byteIndex = position >> 3;
        const uint32_t bitIndex = position & 7;
        const uint32_t bitsToWrite = numBits < 8u - bitIndex ? numBits : 8u - bitIndex;
        const uint32_t mask = (1u << bitsToWrite) - 1u;

        data[byteIndex] |= static_cast<uint8_t>((value & mask) << bitIndex);

        value >>= bitsToWrite;
        position += bitsToWrite;
        numBits -= bitsToWrite;
    }
}

YBI_CPU_GPU int FloatAsInt(float value)
{
#if defined(__CUDA_ARCH__)
    return __float_as_int(value);
#else
    int bits;
    memcpy(&bits, &value, sizeof(bits));
    return bits;
#endif
}

YBI_CPU_GPU float IntAsFloat(int value)
{
#if defined(__CUDA_ARCH__)
    return __int_as_float(value);
#else
    float bits;
    memcpy(&bits, &value, sizeof(value));
    return bits;
#endif
}

YBI_CPU_GPU uint32_t FloatAsUint(float value)
{
#if defined(__CUDA_ARCH__)
    return __float_as_uint(value);
#else
    uint32_t bits;
    memcpy(&bits, &value, sizeof(bits));
    return bits;
#endif
}

YBI_CPU_GPU float UintAsFloat(uint32_t value)
{
#if defined(__CUDA_ARCH__)
    return __uint_as_float(value);
#else
    float bits;
    memcpy(&bits, &value, sizeof(value));
    return bits;
#endif
}

YBI_CPU_GPU uint32_t FloatToOrderedBits(float value)
{
    uint32_t bits = FloatAsUint(value);
    uint32_t mask = (0u - (bits >> 31)) | 0x80000000u;
    return bits ^ mask;
}

YBI_CPU_GPU uint32_t ZigZagEncode(int32_t value)
{
    const uint32_t bits = static_cast<uint32_t>(value);
    return (bits << 1) ^ (0u - (bits >> 31));
}

YBI_CPU_GPU int32_t ZigZagDecode(uint32_t value)
{
    return static_cast<int32_t>((value >> 1) ^ (0u - (value & 1u)));
}

} // namespace ybi
