#pragma once

#include "core/vec2f.h"
#include "core/vec3f.h"

namespace ybi
{

#define MAX_CLOSURE_SIZE 40

struct ClosureEvaluationInfo
{
    Vec3f shadingNormal;
    Vec3f dpdu;
    float *closureData;
};

struct ClosureInitializationInfo : ClosureEvaluationInfo
{
    Vec2f uv;
    Vec3f rayDirection;
    float rayConeWidth;
    float triangleArea;
    float uvArea;
};

} // namespace ybi
