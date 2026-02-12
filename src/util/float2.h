#pragma once

#include "util/assert.h"

YBI_NAMESPACE_BEGIN

struct float2
{
    float x, y;
    __forceinline float operator[](int i) const
    {
        YBI_ASSERT(i >= 0);
        YBI_ASSERT(i < 2);
        return *(&x + i);
    }
    __forceinline float &operator[](int i)
    {
        YBI_ASSERT(i >= 0);
        YBI_ASSERT(i < 2);
        return *(&x + i);
    }
};

__forceinline float2 make_float2(const float x, const float y)
{
    return {x, y};
}

YBI_NAMESPACE_END
