#pragma once

#include "core/math.h"
#include "device/cuda_compat.h"

namespace ybi
{

// NOTE: these are taken from MaterialX implementations
YBI_CPU_GPU float GGXNDF(const Vec3f &h, const Vec2f &alpha)
{
    const Vec2f &he = h.Xy() / alpha;
    float denom = Dot(he, he) + Sqr(h.z);
    return 1.f / (YBI_PI * alpha.x * alpha.y * Sqr(denom));
}

// Height-correlated Smith masking-shadowing
// http://jcgt.org/published/0003/02/03/paper.pdf
// Equations 72 and 99
YBI_CPU_GPU float GGXSmithG2(float nDotL, float nDotV, float alpha)
{
    float alpha2 = Sqr(alpha);
    float lambdaL = Sqrt(alpha2 + (1.f - alpha2) * Sqr(nDotL));
    float lambdaV = Sqrt(alpha2 + (1.f - alpha2) * Sqr(nDotV));
    return 2.f / (lambdaL / nDotL + lambdaV / nDotV);
}

// MaterialX mx_ggx_smith_G1 from mx_microfacet_specular.glsl.
YBI_CPU_GPU float GGXSmithG1(float cosTheta, float alpha)
{
    const float cosTheta2 = Sqr(cosTheta);
    const float tanTheta2 = (1.f - cosTheta2) / cosTheta2;
    return 2.f / (1.f + Sqrt(1.f + Sqr(alpha) * tanTheta2));
}

YBI_CPU_GPU Vec3f FresnelConductor(float cosTheta, Vec3f ior, Vec3f extinction)
{
    float cos2Theta = Sqr(cosTheta);
    float sin2Theta = Max(0.f, 1.f - cos2Theta);

    Vec3f eta2 = Sqr(ior);
    Vec3f etak2 = Sqr(extinction);
    Vec3f intermediateValue0 = eta2 - etak2 - sin2Theta;

    Vec3f a2plusb2 = Sqrt(Sqr(intermediateValue0) + 4.f * eta2 * etak2);
    Vec3f a = Sqrt(0.5f * (a2plusb2 + intermediateValue0));

    Vec3f intermediateValue1 = 2.f * a * cosTheta;
    Vec3f Rs =
        (a2plusb2 + cos2Theta - intermediateValue1) / (a2plusb2 + cos2Theta + intermediateValue1);

    Vec3f intermediateValue2 = cos2Theta * a2plusb2 + Sqr(sin2Theta);
    Vec3f intermediateValue3 = intermediateValue1 * sin2Theta;
    Vec3f Rp =
        Rs * (intermediateValue2 - intermediateValue3) / (intermediateValue2 + intermediateValue3);

    return 0.5f * (Rs + Rp);
}

} // namespace ybi
