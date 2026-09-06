#pragma once

#include "render/bsdf/openpbr/ggx.h"
#include "render/bsdf/openpbr/glossy_diffuse.h"

namespace ybi
{

struct Aggregate
{
    Vec3f baseColor;
    Vec3f baseDiffuseRoughness;

    float specularRoughness;
};

}
