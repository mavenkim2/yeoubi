#pragma once

#include "core/assert.h"
#include "device/cuda_compat.h"

namespace ybi
{

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

} // namespace ybi
