#pragma once

#include "core/vec3f.h"
#include "render/bsdf/materialx/material_defines.h"

namespace ybi
{

#define SHEEN_TINT_OFFSET 0
#define SHEEN_ROUGHNESS_OFFSET (SHEEN_TINT_OFFSET + MATERIAL_NUM_CLOSURE_BITS)
#define SHEEN_PROB_INDEX_OFFSET (SHEEN_ROUGHNESS_OFFSET + MATERIAL_NUM_CLOSURE_BITS)

struct SheenBSDF
{
    Vec3f tint;
    float roughness;

    int probabilityIndex;
};

} // namespace ybi
