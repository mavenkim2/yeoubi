#pragma once

#include "util/assert.h"

YBI_NAMESPACE_BEGIN

struct float3
{
    float x, y, z;
    __forceinline float operator[](int i) const
    {
        YBI_ASSERT(i >= 0);
        YBI_ASSERT(i < 3);
        return *(&x + i);
    }
    __forceinline float &operator[](int i)
    {
        YBI_ASSERT(i >= 0);
        YBI_ASSERT(i < 3);
        return *(&x + i);
    }
};

__forceinline float3 make_float3(const float x, const float y, const float z)
{
    return {x, y, z};
}

YBI_NAMESPACE_END
