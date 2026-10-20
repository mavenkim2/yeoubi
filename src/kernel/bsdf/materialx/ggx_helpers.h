#pragma once

#include "core/math.h"
#include "core/vec2f.h"
#include "core/vec3f.h"
#include "core/vec4f.h"
#include "device/cuda_compat.h"

namespace ybi
{

YBI_CPU_GPU bool IsSmoothMicrofacet(const Vec2f &roughness)
{
    return roughness.x == 0.f && roughness.y == 0.f;
}

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

YBI_CPU_GPU Vec3f GGXDirectionalAlbedo(float nDotV,
                                      float alpha,
                                      const Vec3f &f0,
                                      const Vec3f &f90)
{
    const float x = nDotV;
    const float y = alpha;
    const float x2 = Sqr(x);
    const float y2 = Sqr(y);
    const Vec4f r =
        Vec4f(0.1003, 0.9345, 1.0, 1.0) + Vec4f(-0.6303, -2.323, -1.765, 0.2281) * x +
        Vec4f(9.748, 2.229, 8.263, 15.94) * y + Vec4f(-2.038, -3.748, 11.53, -55.83) * x * y +
        Vec4f(29.34, 1.424, 28.96, 13.08) * x2 + Vec4f(-8.245, -0.7684, -7.507, 41.26) * y2 +
        Vec4f(-26.44, 1.436, -36.11, 54.9) * x2 * y + Vec4f(19.99, 0.2913, 15.86, 300.2) * x * y2 +
        Vec4f(-5.448, 0.6286, 33.37, -285.1) * x2 * y2;
    const Vec2f ab = Clamp(r.Xy() / r.Zw(), 0.f, 1.f);
    return f0 * ab.x + f90 * ab.y;
}

// https://blog.selfshadow.com/publications/turquin/ms_comp_final.pdf
// Equations 14 and 16.
YBI_CPU_GPU Vec3f GGXEnergyCompensation(float nDotV, float alpha, const Vec3f &fss)
{
    const float ess = GGXDirectionalAlbedo(nDotV, alpha, Vec3f(1.f), Vec3f(1.f)).x;
    return 1.f + fss * ((1.f - ess) / Max(ess, 1e-8f));
}

YBI_CPU_GPU Vec3f FresnelConductor(float cosTheta, Vec3f ior, Vec3f extinction)
{
    float cos2Theta = Sqr(cosTheta);
    float sin2Theta = Max(0.f, 1.f - cos2Theta);

    Vec3f eta2 = Sqr(ior);
    Vec3f etak2 = Sqr(extinction);
    Vec3f intermediateValue0 = eta2 - etak2 - sin2Theta;

    Vec3f a2plusb2 = Sqrt(Sqr(intermediateValue0) + 4.f * eta2 * etak2);
    Vec3f a = Sqrt(Max(0.5f * (a2plusb2 + intermediateValue0), 0.f));

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
