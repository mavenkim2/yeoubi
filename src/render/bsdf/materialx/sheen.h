#pragma once

#include "core/math.h"
#include "device/cuda_compat.h"
#include "render/bsdf/materialx/bsdf_types.h"
#include "render/bsdf/materialx/closure.h"
#include "scene/device_scene.h"

namespace ybi
{

YBI_CPU_GPU bool GetSheenParameters(const MaterialDescriptor &descriptor,
                                    const BsdfClosure &closure,
                                    SheenBsdf &bsdf)
{
    const uint32_t base = descriptor.parameterOffset;
    if (!closure.CanLoad(base, SheenParameterCount))
        return false;
    return closure.LoadFloat(base + SheenWeight, bsdf.weight) &&
           closure.LoadFloat3(base + SheenColor, bsdf.color) &&
           closure.LoadFloat(base + SheenRoughness, bsdf.roughness) &&
           closure.LoadFloat3(base + SheenNormal, bsdf.normal);
}

// NOTE: taken from materialx implementation
YBI_CPU_GPU float ImageWorksSheenDirectionalAlbedoAnalytic(float nDotV, float roughness)
{
    Vec2f r = Vec2f(13.67300, 1.0) + Vec2f(-68.78018, 61.57746) * nDotV +
              Vec2f(799.08825, 442.78211) * roughness +
              Vec2f(-905.00061, 2597.49308) * nDotV * roughness +
              Vec2f(60.28956, 121.81241) * Sqr(nDotV) +
              Vec2f(1086.96473, 3045.55075) * Sqr(roughness);
    return r.x / r.y;
}

YBI_CPU_GPU float SheenAlbedo(float roughness, float nDotV)
{
    return Clamp(ImageWorksSheenDirectionalAlbedoAnalytic(nDotV, roughness), 0.f, 1.f);
}

YBI_CPU_GPU Vec3f SampleCharlieHalfVector(Vec2f u, float roughness)
{
    const float phi = 2.0f * YBI_PI * u.x;
    const float sinTheta = Pow(u.y, roughness / (1.f + 2.f * roughness));
    const float cosTheta = Sqrt(Max(0.f, 1.f - sinTheta * sinTheta));

    return Vec3f(sinTheta * Cos(phi), sinTheta * Sin(phi), cosTheta);
}

YBI_CPU_GPU float CharlieNdf(float nDotH, float roughness)
{
    const float invRoughness = 1.f / roughness;
    const float sin2 = Max(0.f, 1.f - nDotH * nDotH);

    return (2.f + invRoughness) * Pow(sin2, 0.5f * invRoughness) / (2.f * YBI_PI);
}

YBI_CPU_GPU bool SampleContyKullaSheen(const Vec3f &viewDirection,
                                       const Vec2f &u,
                                       float roughness,
                                       Vec3f &lightDirection)
{
    const Vec3f h = SampleCharlieHalfVector(u, roughness);

    float hDotV = Dot(h, viewDirection);

    if (hDotV <= 0.f)
    {
        return false;
    }

    lightDirection = -viewDirection + 2.f * hDotV * h;

    float hDotL = Dot(h, lightDirection);

    if (hDotL <= 0.f)
    {
        return false;
    }

#if 0
    float nDotL = lightDirection.z;
    float hDotV = Dot(h, viewDirection);

    float denom = nDotL + nDotV - nDotL * nDotV;
    float mcEstimate = nDotL * hDotV / (h.z * denom);
    weight = mcEstimate * sheenWeight * sheenTint;
#endif

    return true;
}

YBI_CPU_GPU Vec3f EvaluateContyKullaSheen(const Vec3f &viewDirection,
                                          const Vec3f &lightDirection,
                                          const Vec3f &sheenTint,
                                          float sheenWeight,
                                          float roughness)
{
    if (viewDirection.z * lightDirection.z <= 0.f)
    {
        return Vec3f(0.f);
    }

    Vec3f h = Normalize(viewDirection + lightDirection);
    float D = CharlieNdf(h.z, roughness);

    float nDotV = viewDirection.z;
    float nDotL = lightDirection.z;
    float denom = 4.f * (nDotL + nDotV - nDotL * nDotV);

    Vec3f result = sheenTint * sheenWeight * D * nDotL / denom;
    return result;
}

YBI_CPU_GPU float
PDFContyKullaSheen(const Vec3f &viewDirection, const Vec3f &lightDirection, float roughness)
{
    if (viewDirection.z * lightDirection.z <= 0.f)
    {
        return 0.f;
    }

    Vec3f h = Normalize(viewDirection + lightDirection);
    float D = CharlieNdf(h.z, roughness);
    float pdf = D * h.z / (4.f * Dot(h, viewDirection));
    return pdf;
}

// YBI_CPU_GPU Vec3f SampleSheen(const SheenBsdf &bsdf) {}

// Conty kulla sheen
} // namespace ybi
