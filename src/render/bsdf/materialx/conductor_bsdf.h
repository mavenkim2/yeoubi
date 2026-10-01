#pragma once

#include "core/math.h"
#include "core/vec2f.h"
#include "core/vec3f.h"
#include "render/bsdf/materialx/bsdf_types.h"
#include "render/bsdf/materialx/closure.h"
#include "render/bsdf/materialx/dielectric_bsdf.h"
#include "scene/device_scene.h"
#include <cstdint>

namespace ybi
{

enum ConductorParameterOffset : uint32_t
{
    ConductorWeight = 0,
    ConductorIor = 1,
    ConductorExtinction = 4,
    ConductorRoughness = 7,
    ConductorNormal = 9,
    ConductorParameterCount = 12,
};

YBI_CPU_GPU bool GetConductorParameters(const MaterialDescriptor &descriptor,
                                        const BsdfClosure &closure,
                                        ConductorBsdf &bsdf)
{
    const uint32_t base = descriptor.parameterOffset;
    if (!closure.CanLoad(base, ConductorParameterCount))
        return false;
    return closure.LoadFloat(base + ConductorWeight, bsdf.weight) &&
           closure.LoadFloat3(base + ConductorIor, bsdf.ior) &&
           closure.LoadFloat3(base + ConductorExtinction, bsdf.extinction) &&
           closure.LoadFloat2(base + ConductorRoughness, bsdf.alpha) &&
           closure.LoadFloat3(base + ConductorNormal, bsdf.normal);
}

YBI_CPU_GPU Vec3f ConductorF0(const Vec3f &ior, const Vec3f &extinction)
{
    Vec3f f0;
    for (int channel = 0; channel < 3; ++channel)
    {
        const float n = ior[channel];
        const float k2 = Sqr(extinction[channel]);
        f0[channel] = (Sqr(n - 1.f) + k2) / (Sqr(n + 1.f) + k2);
    }
    return f0;
}

// A directional-albedo estimate for lobe importance, not a conductor BRDF.
YBI_CPU_GPU Vec3f ConductorDirectionalAlbedo(const ConductorBsdf &bsdf, float nDotV)
{
    if (bsdf.weight <= 0.f)
        return Vec3f(0.f);

    const Vec3f f0 = ConductorF0(bsdf.ior, bsdf.extinction);
    const float alpha = Sqrt(Max(1e-8f, bsdf.alpha.x) * Max(1e-8f, bsdf.alpha.y));
    const Vec3f f = FresnelConductor(nDotV, bsdf.ior, bsdf.extinction);
    const Vec3f comp = GGXEnergyCompensation(nDotV, alpha, f);
    return Clamp(
        bsdf.weight * GGXDirectionalAlbedo(nDotV, alpha, f0, Vec3f(1.f)) * comp, 0.f, 1.f);
}

YBI_CPU_GPU bool SampleConductorBsdf(const Vec3f &viewDirection,
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

    if (lightDirection.z <= 0.f)
    {
        return false;
    }

    return true;
}

YBI_CPU_GPU Vec3f EvaluateConductorBsdf(const Vec3f &viewDirection,
                                        const Vec3f &lightDirection,
                                        float weight,
                                        const Vec3f &ior,
                                        const Vec3f &extinction,
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
    const Vec3f f = FresnelConductor(vDotH, ior, extinction);
    const Vec3f compensation = GGXEnergyCompensation(nDotV, averageAlpha, f);

    Vec3f result = d * f * g * compensation * weight / (4.f * nDotV);
    // printf("%f %f %f\n", result.x, result.y, result.z);
    return result;
}

YBI_CPU_GPU float
PDFConductorBsdf(const Vec3f &viewDirection, const Vec3f &lightDirection, const Vec2f &roughness)
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

} // namespace ybi
