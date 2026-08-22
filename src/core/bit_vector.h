#pragma once

#include "core/arena.h"
#include "core/bit.h"

namespace ybi
{

struct BitVector
{
    View<uint32_t> bits;
    uint32_t maxNumBits;
    BitVector() {}
    BitVector(Arena *arena, uint32_t maxNumBits);
    void SetBit(uint32_t bit);
    void UnsetBit(uint32_t bit);
    bool GetBit(uint32_t bit) const;
    void WriteBits(uint32_t &position, uint32_t value, uint32_t numBits);
};

inline BitVector::BitVector(Arena *arena, uint32_t maxNumBits) : maxNumBits(maxNumBits)
{
    bits = arena->PushArrayZero<uint32_t>((maxNumBits + 31) >> 5);
}

inline void BitVector::SetBit(uint32_t bit)
{
    YBI_ERROR(bit < maxNumBits, "bit: %u, maxNumBits: %u\n", bit, maxNumBits);
    bits[bit >> 5] |= 1ull << (bit & 31);
}

inline void BitVector::UnsetBit(uint32_t bit)
{
    YBI_ASSERT(bit < maxNumBits);
    bits[bit >> 5] &= ~(1ull << (bit & 31));
}

inline bool BitVector::GetBit(uint32_t bit) const
{
    YBI_ASSERT(bit < maxNumBits);
    return bits[bit >> 5] & (1 << (bit & 31));
}

inline void BitVector::WriteBits(uint32_t &position, uint32_t value, uint32_t numBits)
{
    YBI_ASSERT(position + numBits <= maxNumBits);
    ybi::WriteBits(bits.data(), position, value, numBits);
}

} // namespace ybi
