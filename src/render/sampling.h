#pragma once

#include "core/math.h"
#include "core/vec2f.h"
#include "device/cuda_compat.h"
#include <cmath>
#include <math.h>

namespace ybi
{

YBI_CPU_GPU Vec2f SampleUniformDisk(Vec2f u)
{
    u = 2.f * u - 1.f;

    if (u.x == 0.f && u.y == 0.f)
        return Vec2f(0.f);

    float r;
    float phi;

    if (Abs(u.x) > Abs(u.y))
    {
        r = u.x;
        phi = YBI_PI_OVER_4 * u.y / u.x;
    }
    else
    {
        r = u.y;
        phi = YBI_PI_OVER_2 - YBI_PI_OVER_4 * (u.x / u.y);
    }

    return Vec2f(r * Cos(phi), r * Sin(phi));
}

// Malley's Method
YBI_CPU_GPU Vec3f SampleCosineHemisphere(Vec2f u)
{
#if 0
    float phi = 2.f * YBI_PI * u.x;
    float cosTheta = sqrtf(u.y);
    float sinTheta = sqrtf(1 - u.y);

    Vec3f result(sinTheta * cosf(phi), sinTheta * sinf(phi), cosTheta);
#else
    Vec2f d = SampleUniformDisk(u);
    float z = sqrtf(Max(0.f, 1 - d.x * d.x - d.y * d.y));
    Vec3f result(d.x, d.y, z);
#endif
    return result;
}

} // namespace ybi
