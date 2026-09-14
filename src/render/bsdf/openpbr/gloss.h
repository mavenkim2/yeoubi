#pragma once

#include "render/bsdf/openpbr/ggx.h"

namespace ybi 
{

YBI_GPU bool SampleGlossLobe(const Vec3f &viewDirection,
                             const Vec2f &a,
                             const Vec2f &u,
                             const float eta,
                             const Vec3f &specularColor,
                             Vec3f &lightDirection,
                             Vec3f &weight)
{
    bool success = SampleGGX(viewDirection, a, u, eta, lightDirection, weight);
    if (success)
    {
        weight *= specularColor;
    }
    return success;
}

YBI_CPU_GPU Vec3f EvaluateGlossLobe(const Vec3f &viewDirection,
                                    const Vec3f &lightDirection,
                                    const Vec3f &specularColor,
                                    float ior,
                                    const Vec2f &a)
{
    if (viewDirection.z * lightDirection.z <= 0.f)
    {
        return Vec3f(0.f);
    }
    return EvaluateGGX(viewDirection, lightDirection, ior, a) * specularColor;
}

}
