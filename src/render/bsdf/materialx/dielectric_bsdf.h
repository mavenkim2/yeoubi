#pragma once

#include "core/math.h"
#include <cmath>

namespace ybi
{

struct DielectricBsdf
{
    Vec2f alpha;
    float weight;
    float ior;
};

YBI_GPU void GetDielectricParameters(const DielectricBsdf &dielectric) {}

// NOTE: these are taken from MaterialX implementations
YBI_GPU float GGXNDF(const Vec3f &h, const Vec2f &alpha)
{
    const Vec2f &he = h.Xy() / alpha;
    float denom = Dot(he, he) + Sqr(h.z);
    return 1.f / (YBI_PI * alpha.x * alpha.y * Sqr(denom));
}

// Height-correlated Smith masking-shadowing
// http://jcgt.org/published/0003/02/03/paper.pdf
// Equations 72 and 99
YBI_GPU float GGXSmithG2(float nDotL, float nDotV, float alpha)
{
    float alpha2 = Sqr(alpha);
    float lambdaL = Sqrt(alpha2 + (1.f - alpha2) * Sqr(nDotL));
    float lambdaV = Sqrt(alpha2 + (1.f - alpha2) * Sqr(nDotV));
    return 2.f / (lambdaL / nDotL + lambdaV / nDotV);
}

// https://seblagarde.wordpress.com/2013/04/29/memo-on-fresnel-equations/
YBI_GPU float FresnelDielectric(float cosTheta, float ior)
{
    float c = cosTheta;
    float g2 = ior * ior + c * c - 1.0;
    if (g2 < 0.0)
    {
        return 1.0;
    }

    float g = sqrtf(g2);
    return 0.5 * Sqr((g - c) / (g + c)) * (1.0 + Sqr(((g + c) * c - 1.f) / ((g - c) * c + 1.f)));
}

YBI_GPU float IorToF0(float ior)
{
    return Sqr((ior - 1.0) / (ior + 1.0));
}

YBI_GPU Vec3f DielectricBsdfEvaluate(const DielectricBsdf &dielectric,
                                     const Vec3f &lightDirection,
                                     const Vec3f &viewDirection,
                                     float ior)
{
    if (lightDirection.z * viewDirection.z < 0.f)
    {
        return Vec3f(0.f);
    }

    Vec3f h = Normalize(lightDirection + viewDirection);
    float averageAlpha = Sqrt(dielectric.alpha.x * dielectric.alpha.y);

    const float D = GGXNDF(h, dielectric.alpha);
    const float G = GGXSmithG2(lightDirection.z, viewDirection.z, averageAlpha);
    const float F = FresnelDielectric(Dot(viewDirection, h), ior);

    const float eval = D * G * F / (4 * viewDirection.z);
    return Vec3f(eval);
}

YBI_GPU Vec3f DielectricBsdfAlbedo(float nDotV, float alpha, const Vec3f &f0, const Vec3f &f90)
{
    float x = nDotV;
    float y = alpha;
    float x2 = Sqr(x);
    float y2 = Sqr(y);
    Vec4f r =
        Vec4f(0.1003, 0.9345, 1.0, 1.0) + Vec4f(-0.6303, -2.323, -1.765, 0.2281) * x +
        Vec4f(9.748, 2.229, 8.263, 15.94) * y + Vec4f(-2.038, -3.748, 11.53, -55.83) * x * y +
        Vec4f(29.34, 1.424, 28.96, 13.08) * x2 + Vec4f(-8.245, -0.7684, -7.507, 41.26) * y2 +
        Vec4f(-26.44, 1.436, -36.11, 54.9) * x2 * y + Vec4f(19.99, 0.2913, 15.86, 300.2) * x * y2 +
        Vec4f(-5.448, 0.6286, 33.37, -285.1) * x2 * y2;
    Vec2f ab = Clamp(r.Xy() / r.Zw(), 0.f, 1.f);
    return f0 * ab.x + f90 * ab.y;
}

// https://blog.selfshadow.com/publications/turquin/ms_comp_final.pdf
// Equations 14 and 16
YBI_GPU Vec3f GGXEnergyCompensation(float nDotV, float alpha, const Vec3f &fss)
{
    Vec3f ess = DielectricBsdfAlbedo(nDotV, alpha, Vec3f(1.f), Vec3f(1.f));
    return 1.0 + fss * (1.0 - ess) / ess;
}

// Compute the cosine-weighted average of the Fresnel reflectance over the hemisphere.
// https://blog.selfshadow.com/publications/s2017-shading-course/imageworks/s2017_pbs_imageworks_slides_v2.pdf
YBI_GPU Vec3f FresnelAverage(float ior)
{
    float f0 = FresnelDielectric(1.f, ior);
    Vec3f f90(1.f);
    return f0 + (f90 - f0) * (1.f / 21.f);
}

YBI_GPU Vec3f CalculateDielectricBsdfAlbedo(const DielectricBsdf &bsdf, const Vec3f &viewDirection)
{
    float alpha = Sqrt(bsdf.alpha.x * bsdf.alpha.y);
    float f0 = IorToF0(bsdf.ior);

    Vec3f fAvg = FresnelAverage(bsdf.ior);
    Vec3f comp = GGXEnergyCompensation(viewDirection.z, alpha, fAvg);
    Vec3f albedo =
        bsdf.weight * DielectricBsdfAlbedo(viewDirection.z, alpha, Vec3f(f0), Vec3f(1.f)) * comp;
    return albedo;
}

} // namespace ybi
