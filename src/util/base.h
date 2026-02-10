#pragma once

#include "util/assert.h"

YBI_NAMESPACE_BEGIN

namespace util
{

__forceinline size_t AlignUp(size_t val, size_t align)
{
    return (val + align - 1) & ~(align - 1);
}

template <typename Dest, typename Src>
__forceinline void Copy(Dest dst, Src src, size_t count)
{
    using DestT = typename Dest::value_type;
    using SrcT = typename Src::value_type;
    static_assert(std::is_same_v<DestT, SrcT>, "Can only copy to and from the same type.\n");

    ASSERT(dst.size() >= count);
    ASSERT(src.size() >= count);

    memcpy(dst.data(), src.data(), sizeof(DestT) * count);
}

} // namespace util

YBI_NAMESPACE_END
