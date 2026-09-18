#pragma once

#include "core/assert.h"
#include "device/cuda_compat.h"
#include "render/bsdf/materialx/material_defines.h"

namespace ybi
{

struct BsdfClosure
{
    float values[MATERIAL_NUM_CLOSURE_VALUES];

    YBI_GPU float LoadFloat(int src)
    {
        YBI_ASSERT(src < MATERIAL_NUM_CLOSURE_VALUES);
        float result = values[src];
        return result;
    }

    YBI_GPU Vec2f LoadFloat2(int src)
    {
        YBI_ASSERT(src + 1 < MATERIAL_NUM_CLOSURE_VALUES);
        Vec2f v;
        v[0] = values[src];
        v[1] = values[src + 1];
        return v;
    }

    YBI_GPU Vec3f LoadFloat3(int src)
    {
        YBI_ASSERT(src + 2 < MATERIAL_NUM_CLOSURE_VALUES);
        Vec3f v;
        v[0] = values[src];
        v[1] = values[src + 1];
        v[2] = values[src + 2];
        return v;
    }

    YBI_GPU Vec4f LoadFloat4(int src)
    {
        YBI_ASSERT(src + 3 < MATERIAL_NUM_CLOSURE_VALUES);
        Vec4f v;
        v[0] = values[src];
        v[1] = values[src + 1];
        v[2] = values[src + 2];
        v[3] = values[src + 3];
        return v;
    }
};

} // namespace ybi
