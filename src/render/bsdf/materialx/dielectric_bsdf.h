#pragma once

#include "core/math.h"
#include "core/vec2f.h"
#include "core/vec3f.h"
#include "core/vec4f.h"
#include "render/bsdf/materialx/closure.h"
#include "scene/device_scene.h"
#include <cmath>
#include <cstdint>

namespace ybi
{

struct DielectricBsdf
{
    Vec2f alpha;
    float weight;
    float ior;
    Vec3f tint;
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

YBI_CPU_GPU bool GetDielectricParameters(const MaterialDescriptor &descriptor,
                                         const BsdfClosure &closure,
                                         DielectricBsdf &bsdf)
{
    const uint32_t base = descriptor.parameterOffset;
    if (!closure.CanLoad(base, DielectricParameterCount))
        return false;
    return closure.LoadFloat(base + DielectricWeight, bsdf.weight) &&
           closure.LoadFloat3(base + DielectricTint, bsdf.tint) &&
           closure.LoadFloat(base + DielectricIor, bsdf.ior) &&
           closure.LoadFloat2(base + DielectricRoughness, bsdf.alpha) &&
           closure.LoadFloat3(base + DielectricNormal, bsdf.normal);
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

// Direct port of MaterialX mx_ggx_importance_sample_VNDF.
YBI_CPU_GPU Vec3f SampleMaterialXGGXVNDF(const Vec2f &u,
                                         const Vec3f &viewDirection,
                                         const Vec2f &alpha)
{
    Vec3f view = Normalize(Vec3f(viewDirection.Xy() * alpha, viewDirection.z));

    const float phi = 2.f * YBI_PI * u.x;
    const float z = (1.f - u.y) * (1.f + view.z) - view.z;
    const float sinTheta = Sqrt(Clamp(1.f - z * z, 0.f, 1.f));
    const Vec3f sphericalCap(sinTheta * Cos(phi), sinTheta * Sin(phi), z);

    Vec3f h = sphericalCap + view;
    h = Normalize(Vec3f(h.Xy() * alpha, Max(h.z, 0.f)));
    return h;
}

// Reflection-only GGX sampling. MaterialX uses this VNDF procedure for its
// filtered-importance sampling path.
YBI_CPU_GPU bool SampleDielectricBsdf(const Vec3f &viewDirection,
                                      const Vec2f &u,
                                      const Vec2f &roughness,
                                      Vec3f &lightDirection)
{
    if (viewDirection.z <= 0.f)
        return false;

    const Vec2f alpha = Clamp(roughness, 1e-8f, 1.f);
    const Vec3f h = SampleMaterialXGGXVNDF(u, viewDirection, alpha);
    const float hDotV = Dot(h, viewDirection);
    if (hDotV <= 0.f)
        return false;

    lightDirection = -viewDirection + 2.f * hDotV * h;
    return lightDirection.z > 0.f;
}

// Reflected-direction PDF used by MaterialX mx_environment_radiance after
// sampling mx_ggx_importance_sample_VNDF.
YBI_CPU_GPU float PDFDielectricBsdf(const Vec3f &viewDirection,
                                    const Vec3f &lightDirection,
                                    const Vec2f &roughness)
{
    if (viewDirection.z <= 0.f || lightDirection.z <= 0.f)
        return 0.f;

    const Vec2f alpha = Clamp(roughness, 1e-8f, 1.f);
    const Vec3f h = Normalize(viewDirection + lightDirection);
    if (Dot(h, viewDirection) <= 0.f)
        return 0.f;

    const float nDotV = Clamp(viewDirection.z, 1e-8f, 1.f);
    const float averageAlpha = Sqrt(alpha.x * alpha.y);
    const float g1V = GGXSmithG1(nDotV, averageAlpha);
    return GGXNDF(h, alpha) * g1V / (4.f * nDotV);
}

// https://seblagarde.wordpress.com/2013/04/29/memo-on-fresnel-equations/
YBI_CPU_GPU float FresnelDielectric(float cosTheta, float ior)
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

YBI_CPU_GPU float IorToF0(float ior)
{
    return Sqr((ior - 1.0) / (ior + 1.0));
}

YBI_CPU_GPU Vec3f DielectricBsdfAlbedo(float nDotV,
                                      float alpha,
                                      const Vec3f &f0,
                                      const Vec3f &f90)
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
YBI_CPU_GPU Vec3f GGXEnergyCompensation(float nDotV, float alpha, const Vec3f &fss)
{
    Vec3f ess = DielectricBsdfAlbedo(nDotV, alpha, Vec3f(1.f), Vec3f(1.f));
    return 1.0 + fss * (1.0 - ess) / Max(ess, 1e-8f);
}

// Compute the cosine-weighted average of the Fresnel reflectance over the hemisphere.
// https://blog.selfshadow.com/publications/s2017-shading-course/imageworks/s2017_pbs_imageworks_slides_v2.pdf
YBI_CPU_GPU Vec3f FresnelAverage(float ior)
{
    float f0 = FresnelDielectric(1.f, ior);
    Vec3f f90(1.f);
    return f0 + (f90 - f0) * (1.f / 21.f);
}

YBI_CPU_GPU Vec3f CalculateDielectricBsdfAlbedo(const DielectricBsdf &bsdf, float nDotV)
{
    if (bsdf.weight <= 0.f)
        return Vec3f(0.f);

    const Vec2f safeAlpha = Clamp(bsdf.alpha, 1e-8f, 1.f);
    nDotV = Clamp(nDotV, 1e-8f, 1.f);
    float alpha = Sqrt(safeAlpha.x * safeAlpha.y);
    float f0 = IorToF0(bsdf.ior);

    Vec3f fAvg = FresnelAverage(bsdf.ior);
    Vec3f comp = GGXEnergyCompensation(nDotV, alpha, fAvg);
    Vec3f albedo =
        bsdf.weight * DielectricBsdfAlbedo(nDotV, alpha, Vec3f(f0), Vec3f(1.f)) * comp;
    return Clamp(albedo, 0.f, 1.f);
}

// MaterialX reflection response. As in mx_dielectric_bsdf.glsl, this returns
// BRDF * NdotL and includes energy compensation, tint, and node weight.
YBI_CPU_GPU Vec3f EvaluateDielectricBsdf(const Vec3f &viewDirection,
                                         const Vec3f &lightDirection,
                                         float weight,
                                         const Vec3f &tint,
                                         float ior,
                                         const Vec2f &roughness)
{
    if (weight < 1e-8f || viewDirection.z <= 0.f || lightDirection.z <= 0.f)
        return Vec3f(0.f);

    const Vec2f alpha = Clamp(roughness, 1e-8f, 1.f);
    const float nDotV = Clamp(viewDirection.z, 1e-8f, 1.f);
    const float nDotL = Clamp(lightDirection.z, 1e-8f, 1.f);
    const float averageAlpha = Sqrt(alpha.x * alpha.y);
    const Vec3f h = Normalize(lightDirection + viewDirection);
    const float vDotH = Clamp(Dot(viewDirection, h), 1e-8f, 1.f);

    const float d = GGXNDF(h, alpha);
    const float g = GGXSmithG2(nDotL, nDotV, averageAlpha);
    const float f = FresnelDielectric(vDotH, ior);
    const Vec3f compensation =
        GGXEnergyCompensation(nDotV, averageAlpha, Vec3f(f));

    return d * f * g * compensation * Max(tint, 0.f) * weight / (4.f * nDotV);
}

YBI_CPU_GPU float DielectricAlbedo(float nDotV,
                                   float weight,
                                   float ior,
                                   const Vec2f &roughness)
{
    DielectricBsdf bsdf;
    bsdf.alpha = roughness;
    bsdf.weight = weight;
    bsdf.ior = ior;
    return CalculateDielectricBsdfAlbedo(bsdf, nDotV).x;
}

} // namespace ybi
