#pragma once

#include "core/vec3f.h"
#include <cstdint>

namespace ybi
{

enum SheenParameterOffset : uint32_t
{
    SheenWeight = 0,
    SheenColor = 1,
    SheenRoughness = 4,
    SheenNormal = 5,
    SheenParameterCount = 8,
};

struct SheenBSDF
{
    float weight;
    Vec3f color;
    float roughness;
    Vec3f normal;
};

} // namespace ybi
