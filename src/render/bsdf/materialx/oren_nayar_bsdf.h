#pragma once

#include "core/math.h"
#include "core/vec2f.h"
#include "core/vec3f.h"
#include "render/bsdf/materialx/closure.h"
#include "render/sampling.h"
#include "scene/device_scene.h"
#include <cstdint>

namespace ybi
{

enum OrenNayarParameterOffset : uint32_t
{
    OrenNayarWeight = 0,
    OrenNayarColor = 1,
    OrenNayarRoughness = 4,
    OrenNayarNormal = 5,
    OrenNayarParameterCount = 8,
};

struct OrenNayarBsdf
{
    float weight;
    Vec3f color;
    float roughness;
    Vec3f normal;
};

YBI_CPU_GPU bool GetOrenNayarParameters(const MaterialDescriptor &descriptor,
                                        const BsdfClosure &closure,
                                        OrenNayarBsdf &bsdf)
{
    const uint32_t base = descriptor.parameterOffset;
    if (!closure.CanLoad(base, OrenNayarParameterCount))
        return false;
    return closure.LoadFloat(base + OrenNayarWeight, bsdf.weight) &&
           closure.LoadFloat3(base + OrenNayarColor, bsdf.color) &&
           closure.LoadFloat(base + OrenNayarRoughness, bsdf.roughness) &&
           closure.LoadFloat3(base + OrenNayarNormal, bsdf.normal);
}

// MaterialX's qualitative Oren-Nayar model. These functions return BRDF * NdotL,
// matching the other generated MaterialX BSDF evaluators.
YBI_CPU_GPU float EvaluateOrenNayarReflectance(const Vec3f &viewDirection,
                                               const Vec3f &lightDirection,
                                               float roughness)
{
    if (viewDirection.z <= 0.f || lightDirection.z <= 0.f)
        return 0.f;

    const float nDotV = Clamp(viewDirection.z, 1e-8f, 1.f);
    const float nDotL = Clamp(lightDirection.z, 1e-8f, 1.f);
    const float lDotV = Clamp(Dot(lightDirection, viewDirection), 1e-8f, 1.f);
    const float s = lDotV - nDotL * nDotV;
    const float stinv = s > 0.f ? s / Max(nDotL, nDotV) : 0.f;

    const float sigma2 = Sqr(roughness);
    const float a = 1.f - 0.5f * sigma2 / (sigma2 + 0.33f);
    const float b = 0.45f * sigma2 / (sigma2 + 0.09f);
    return a + b * stinv;
}

YBI_CPU_GPU bool
SampleOrenNayarBsdf(const Vec3f &viewDirection, const Vec2f &u, Vec3f &lightDirection)
{
    if (viewDirection.z <= 0.f)
        return false;

    lightDirection = SampleCosineHemisphere(u);
    return true;
}

YBI_CPU_GPU Vec3f EvaluateOrenNayarBsdf(const Vec3f &viewDirection,
                                        const Vec3f &lightDirection,
                                        float weight,
                                        const Vec3f &color,
                                        float roughness)
{
    if (weight < 1e-8f)
        return Vec3f(0.f);

    const float reflectance =
        EvaluateOrenNayarReflectance(viewDirection, lightDirection, roughness);
    return weight * Max(color, 0.f) * reflectance * lightDirection.z * YBI_RCP_PI;
}

YBI_CPU_GPU float PDFOrenNayarBsdf(const Vec3f &viewDirection, const Vec3f &lightDirection)
{
    if (viewDirection.z <= 0.f || lightDirection.z <= 0.f)
        return 0.f;
    return lightDirection.z * YBI_RCP_PI;
}

// MaterialX's analytic directional-albedo fit for uncompensated Oren-Nayar.
YBI_CPU_GPU float OrenNayarDirectionalAlbedo(float nDotV, float roughness)
{
    const float r0 =
        1.f - 0.4297f * roughness - 0.7632f * nDotV * roughness + 1.4385f * Sqr(roughness);
    const float r1 =
        1.f - 0.6076f * roughness - 0.4993f * nDotV * roughness + 2.0315f * Sqr(roughness);
    return Clamp(r0 / r1, 0.f, 1.f);
}

} // namespace ybi
