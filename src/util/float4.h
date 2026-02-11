#pragma once

#include "util/assert.h"

YBI_NAMESPACE_BEGIN

struct float4
{
    float x, y, z, w;
    __forceinline float operator[](int i) const
    {
        YBI_ASSERT(i >= 0);
        YBI_ASSERT(i < 4);
        return *(&x + i);
    }
    __forceinline float &operator[](int i)
    {
        YBI_ASSERT(i >= 0);
        YBI_ASSERT(i < 4);
        return *(&x + i);
    }
};

__forceinline float4 make_float4(const float x, const float y, const float z, const float w)
{
    return {x, y, z, w};
}

YBI_NAMESPACE_END
