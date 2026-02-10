#pragma once

YBI_NAMESPACE_BEGIN

namespace util
{

__forceinline size_t AlignUp(size_t val, size_t align)
{
    return (val + align - 1) & ~(align - 1);
}

} // namespace util

YBI_NAMESPACE_END
