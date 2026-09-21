#pragma once

#include "core/vec2f.h"
#include "core/vec3f.h"
#include "core/vec4f.h"
#include "device/cuda_compat.h"
#include <cstdint>

namespace ybi
{

struct BsdfClosure
{
    const float *values = nullptr;
    uint32_t valueCount = 0;

    YBI_CPU_GPU bool CanLoad(uint32_t src, uint32_t count) const
    {
        return values != nullptr && src <= valueCount && count <= valueCount - src;
    }

    YBI_CPU_GPU bool LoadFloat(uint32_t src, float &result) const
    {
        if (!CanLoad(src, 1))
            return false;
        result = values[src];
        return true;
    }

    YBI_CPU_GPU bool LoadFloat2(uint32_t src, Vec2f &result) const
    {
        if (!CanLoad(src, 2))
            return false;
        result = Vec2f(values[src], values[src + 1]);
        return true;
    }

    YBI_CPU_GPU bool LoadFloat3(uint32_t src, Vec3f &result) const
    {
        if (!CanLoad(src, 3))
            return false;
        result = Vec3f(values[src], values[src + 1], values[src + 2]);
        return true;
    }

    YBI_CPU_GPU bool LoadFloat4(uint32_t src, Vec4f &result) const
    {
        if (!CanLoad(src, 4))
            return false;
        result = Vec4f(values[src], values[src + 1], values[src + 2], values[src + 3]);
        return true;
    }
};

} // namespace ybi
