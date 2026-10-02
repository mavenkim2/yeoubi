#pragma once

#include "core/math.h"
#include "core/vec2f.h"
#include "core/vec3f.h"
#include "device/cuda_compat.h"
#include "kernel/bsdf/materialx/bsdf_types.h"
#include "kernel/bsdf/materialx/closure.h"
#include "kernel/bsdf/materialx/ggx_helpers.h"
#include "kernel/scene.h"

namespace ybi
{

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
YBI_CPU_GPU float
PDFDielectricBsdf(const Vec3f &viewDirection, const Vec3f &lightDirection, const Vec2f &roughness)
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

YBI_CPU_GPU Vec3f CalculateDielectricBsdfAlbedo(const DielectricBsdf &bsdf, float nDotV)
{
    if (bsdf.weight <= 0.f)
        return Vec3f(0.f);

    const Vec2f safeAlpha = Clamp(bsdf.alpha, 1e-8f, 1.f);
    nDotV = Clamp(nDotV, 1e-8f, 1.f);
    float alpha = Sqrt(safeAlpha.x * safeAlpha.y);
    float f0 = IorToF0(bsdf.ior);

    Vec3f f = Vec3f(FresnelDielectric(nDotV, bsdf.ior));
    Vec3f comp = GGXEnergyCompensation(nDotV, alpha, f);
    Vec3f albedo = bsdf.weight * GGXDirectionalAlbedo(nDotV, alpha, Vec3f(f0), Vec3f(1.f)) * comp;
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
    const Vec3f compensation = GGXEnergyCompensation(nDotV, averageAlpha, Vec3f(f));

    return d * f * g * compensation * Max(tint, 0.f) * weight / (4.f * nDotV);
}

YBI_CPU_GPU float DielectricAlbedo(float nDotV, float weight, float ior, const Vec2f &roughness)
{
    DielectricBsdf bsdf;
    bsdf.alpha = roughness;
    bsdf.weight = weight;
    bsdf.ior = ior;
    return CalculateDielectricBsdfAlbedo(bsdf, nDotV).x;
}

} // namespace ybi
