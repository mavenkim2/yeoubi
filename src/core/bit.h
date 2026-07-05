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

} // namespace ybi
