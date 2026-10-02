#pragma once

#include "core/vec2f.h"
#include "core/vec3f.h"
#include "kernel/bsdf/closure_types.h"

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

struct SheenBsdf
{
    float weight;
    Vec3f color;
    float roughness;
    Vec3f normal;
};

struct ConductorBsdf
{
    float weight;
    Vec3f ior;
    Vec3f extinction;
    Vec2f alpha;
    Vec3f normal;
};

enum DielectricParameterOffset : uint32_t
{
    DielectricWeight = 0,
    DielectricTint = 1,
    DielectricIor = 4,
    DielectricRoughness = 5,
    DielectricNormal = 7,
    DielectricParameterCount = 10,
};

struct DielectricBsdf
{
    Vec2f alpha;
    float weight;
    float ior;
    Vec3f tint;
    Vec3f normal;
};

enum OrenNayarParameterOffset : uint32_t
{
    OrenNayarWeight = 0,
    OrenNayarColor = 1,
    OrenNayarRoughness = 4,
    OrenNayarNormal = 5,
    OrenNayarParameterCount = 8,
};

struct OrenNayarBsdf
{
    float weight;
    Vec3f color;
    float roughness;
    Vec3f normal;
};

} // namespace ybi
