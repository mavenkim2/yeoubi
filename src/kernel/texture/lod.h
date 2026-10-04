#pragma once

#include "core/math.h"
#include "core/vec3f.h"
#include "device/cuda_compat.h"

namespace ybi
{

YBI_GPU float ComputeTextureLOD(float triangleArea,
                                float textureArea,
                                Vec3f dir,
                                Vec3f normal,
                                float width,
                                int textureWidth,
                                int textureHeight)
{
    const float safeTextureArea = Max(textureArea, 1.0e-20f);
    const float safeTriangleArea = Max(triangleArea, 1.0e-20f);
    const float safeWidth = Max(Abs(width), 1.0e-6f);
    const float safeNormalDot = Max(Abs(Dot(dir, normal)), 1.0e-4f);

    float lambda = 0.f;
    lambda += 0.5f * log2f(safeTextureArea / safeTriangleArea);
    lambda += log2f(safeWidth);
    lambda += 0.5f * log2f(float(textureWidth) * float(textureHeight));
    lambda -= log2f(safeNormalDot);

    return lambda;
}

} // namespace ybi
