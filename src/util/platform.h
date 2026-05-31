#pragma once

namespace ybi
{

inline uint32_t FirstBitLow(uint32_t x)
{
    if (x == 0)
        return ~0u;
#ifdef _WIN32

#if defined(__GNUC__) || defined(__clang__)
    return (uint32_t)__builtin_ctz(value);
#elif defined(_MSC_VER)
    unsigned long index;
    _BitScanForward(&index, value);
    return (uint32_t)index;
#else
    uint32_t count = 0;
    if ((value & 0xFFFF) == 0)
    {
        value >>= 16;
        count += 16;
    }
    if ((value & 0xFF) == 0)
    {
        value >>= 8;
        count += 8;
    }
    if ((value & 0xF) == 0)
    {
        value >>= 4;
        count += 4;
    }
    if ((value & 0x3) == 0)
    {
        value >>= 2;
        count += 2;
    }
    if ((value & 0x1) == 0)
    {
        count += 1;
    }
    return count;
#endif
#endif
}

} // namespace ybi
