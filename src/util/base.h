#pragma once

#include "util/assert.h"
#include "util/forceinline.h"

namespace ybi
{

namespace util
{

__forceinline size_t AlignUp(size_t val, size_t align)
{
    return (val + align - 1) & ~(align - 1);
}

} // namespace util

} // namespace ybi
