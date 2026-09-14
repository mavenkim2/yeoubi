#pragma once

#include "core/frame.h"
#include "core/vec2f.h"
#include "core/vec3f.h"

namespace ybi
{

enum LobeType
{
    THIN_WALL_TRANS,
    DIFFUSE,
    SPECULAR,
    COUNT,
};

struct OpenPbrBsdf
{
    OrthonormalBasis tangentFrame{Vec3f(1.f, 0.f, 0.f),
                                  Vec3f(0.f, 1.f, 0.f),
                                  Vec3f(0.f, 0.f, 1.f)};

    float fuzzWeight;
    Vec3f fuzzColor;
    float fuzzRoughness;

    float specularReflectionIor;
    float specularRefractionIor;
    Vec2f specularRoughness;
    Vec3f specularColor;

    Vec3f baseColor;

    // Thin wall
    Vec3f thinWallTransmissionColor;
    bool thinWall;

    bool inside;

    // Coat
    Vec3f coatColor;
    float coatWeight;
    Vec2f coatRoughness;
    Vec3f coatDarkeningFactor;
    float coatEnergyReflected;
    float coatIor;

    float lobeWeights[COUNT];
};

} // namespace ybi
