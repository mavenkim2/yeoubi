#pragma once

#include "render/integrator_common.h"
#include "scene/material_light.h"

namespace ybi
{
namespace render
{
namespace integrator
{

struct EvaluatedMaterial
{
    Vec3 baseColor = Vec3(0.7f, 0.7f, 0.7f);
    Vec3 emissiveColor = Vec3(0.0f, 0.0f, 0.0f);
    float roughness = 1.0f;
    float metallic = 0.0f;
    float ior = 1.5f;
    float opacity = 1.0f;
    unsigned int flags = 0u;
};

struct BsdfSample
{
    Vec3 wi = Vec3(0.0f, 0.0f, 1.0f);
    Vec3 f = Vec3(0.0f, 0.0f, 0.0f);
    float pdf = 0.0f;
    bool skipNeeMis = false;
};

YBI_INTEGRATOR_HD float MaxFloat(float a, float b)
{
    return a > b ? a : b;
}

YBI_INTEGRATOR_HD float MinFloat(float a, float b)
{
    return a < b ? a : b;
}

YBI_INTEGRATOR_HD Vec3 MulComponents(const Vec3 &a, const Vec3 &b)
{
    return a * b;
}

YBI_INTEGRATOR_HD Vec3 ClampVec3(const Vec3 &v, float lo, float hi)
{
    return Clamp(v, lo, hi);
}

YBI_INTEGRATOR_HD float Saturate(float value)
{
    return Clamp01(value);
}

YBI_INTEGRATOR_HD float Pow5(float x)
{
    const float x2 = x * x;
    return x2 * x2 * x;
}

YBI_INTEGRATOR_HD float SafeSqrt(float x)
{
    return sqrtf(MaxFloat(x, 0.0f));
}

YBI_INTEGRATOR_HD float ComputeAlpha(const EvaluatedMaterial &material)
{
    return MaxFloat(material.roughness * material.roughness, 1.0e-6f);
}

YBI_INTEGRATOR_HD bool ShouldSkipNeeMis(const EvaluatedMaterial &material)
{
    return material.roughness <= 1.0e-3f;
}

YBI_INTEGRATOR_HD Vec3 FresnelSchlick(float cosTheta, const Vec3 &f0)
{
    const float m = Pow5(1.0f - Saturate(cosTheta));
    return Add(f0, Mul(Sub(Vec3(1.0f, 1.0f, 1.0f), f0), m));
}

YBI_INTEGRATOR_HD float Dggx(float alpha, float nDotH)
{
    const float a2 = alpha * alpha;
    const float denom = MaxFloat((nDotH * nDotH) * (a2 - 1.0f) + 1.0f, 1.0e-6f);
    return a2 / (3.14159265358979323846f * denom * denom);
}

YBI_INTEGRATOR_HD float GSchlickGgx(float alpha, float nDotV)
{
    const float k = 0.5f * alpha;
    return nDotV / MaxFloat(nDotV * (1.0f - k) + k, 1.0e-6f);
}

YBI_INTEGRATOR_HD float GSmith(float alpha, float nDotV, float nDotL)
{
    return GSchlickGgx(alpha, nDotV) * GSchlickGgx(alpha, nDotL);
}

YBI_INTEGRATOR_HD Vec3 SampleGgxHalfVector(float alpha, float u1, float u2)
{
    const float a2 = alpha * alpha;
    const float phi = 6.28318530718f * u2;
    const float tan2Theta = a2 * u1 / MaxFloat(1.0f - u1, 1.0e-6f);
    const float cosTheta = 1.0f / SafeSqrt(1.0f + tan2Theta);
    const float sinTheta = SafeSqrt(1.0f - cosTheta * cosTheta);
    return Vec3(cosf(phi) * sinTheta, sinf(phi) * sinTheta, cosTheta);
}

YBI_INTEGRATOR_HD Vec3 LocalToWorld(const Vec3 &local, const Vec3 &t, const Vec3 &b, const Vec3 &n)
{
    return Add(Add(Mul(t, local.x), Mul(b, local.y)), Mul(n, local.z));
}

YBI_INTEGRATOR_HD Vec3 ComputeF0(const EvaluatedMaterial &material)
{
    const float dielectric = (material.ior - 1.0f) / MaxFloat(material.ior + 1.0f, 1.0e-6f);
    const float dielectricF0 = dielectric * dielectric;
    return Lerp(
        Vec3(dielectricF0, dielectricF0, dielectricF0), material.baseColor, material.metallic);
}

YBI_INTEGRATOR_HD float ComputeSpecularProbability(const EvaluatedMaterial &material,
                                                   const Vec3 &f0)
{
    const float diffuseWeight =
        (1.0f - material.metallic) * MaxFloat(Luminance(material.baseColor), 0.0f);
    const float specularWeight = 0.02f + MaxFloat(Luminance(f0), 0.0f);
    if (diffuseWeight <= 1.0e-6f)
    {
        return 1.0f;
    }
    return Clamp01(specularWeight / MaxFloat(diffuseWeight + specularWeight, 1.0e-6f));
}

YBI_INTEGRATOR_HD float ComputeDiffusePdf(float nDotL)
{
    return nDotL * 0.31830988618f;
}

YBI_INTEGRATOR_HD float ComputeSpecularPdf(const EvaluatedMaterial &material,
                                           const Vec3 &normal,
                                           const Vec3 &wo,
                                           const Vec3 &wi)
{
    const Vec3 halfUnnormalized = Add(wo, wi);
    const float halfLen2 = Dot(halfUnnormalized, halfUnnormalized);
    if (halfLen2 <= 1.0e-8f)
    {
        return 0.0f;
    }

    const Vec3 halfVec = Normalize(halfUnnormalized);
    const float nDotH = MaxFloat(Dot(normal, halfVec), 0.0f);
    const float vDotH = MaxFloat(Dot(wo, halfVec), 0.0f);
    if (nDotH <= 0.0f || vDotH <= 0.0f)
    {
        return 0.0f;
    }

    const float d = Dggx(ComputeAlpha(material), nDotH);
    return d * nDotH / MaxFloat(4.0f * vDotH, 1.0e-6f);
}

YBI_INTEGRATOR_HD Vec3 EvaluateBsdf(const EvaluatedMaterial &material,
                                    const Vec3 &normal,
                                    const Vec3 &wo,
                                    const Vec3 &wi,
                                    float *outPdf = nullptr)
{
    const float nDotV = MaxFloat(Dot(normal, wo), 0.0f);
    const float nDotL = MaxFloat(Dot(normal, wi), 0.0f);
    if (nDotV <= 0.0f || nDotL <= 0.0f)
    {
        if (outPdf)
        {
            *outPdf = 0.0f;
        }
        return Vec3(0.0f, 0.0f, 0.0f);
    }

    const Vec3 halfUnnormalized = Add(wo, wi);
    if (Dot(halfUnnormalized, halfUnnormalized) <= 1.0e-8f)
    {
        if (outPdf)
        {
            *outPdf = 0.0f;
        }
        return Vec3(0.0f, 0.0f, 0.0f);
    }

    const Vec3 halfVec = Normalize(halfUnnormalized);
    const float nDotH = MaxFloat(Dot(normal, halfVec), 0.0f);
    const float vDotH = MaxFloat(Dot(wo, halfVec), 0.0f);

    const Vec3 f0 = ComputeF0(material);
    const Vec3 fresnel = FresnelSchlick(vDotH, f0);
    const float alpha = ComputeAlpha(material);
    const float d = Dggx(alpha, nDotH);
    const float g = GSmith(alpha, nDotV, nDotL);
    const Vec3 specular = Mul(fresnel, d * g / MaxFloat(4.0f * nDotV * nDotL, 1.0e-6f));

    const Vec3 diffuseColor = Mul(material.baseColor, 1.0f - material.metallic);
    const Vec3 diffuse = Mul(diffuseColor, 0.31830988618f);

    if (outPdf)
    {
        const float specProb = ComputeSpecularProbability(material, f0);
        const float diffuseProb = 1.0f - specProb;
        const float diffusePdf = ComputeDiffusePdf(nDotL);
        const float specularPdf = ComputeSpecularPdf(material, normal, wo, wi);
        *outPdf = diffuseProb * diffusePdf + specProb * specularPdf;
    }

    return Add(diffuse, specular);
}

// YBI_INTEGRATOR_HD void EvaluateOrenNayar(float3 wi, float3 wo)
// {
//     float s = ωi⋅ωo−(n⋅ωi)(n⋅ωo) = cos(ϕi−ϕo) sinθisinθo;
// }

YBI_INTEGRATOR_HD void SampleOrenNayar(float3 help)
{
    float3 h = Normalize(help);
}

YBI_INTEGRATOR_HD bool SampleBsdf(const EvaluatedMaterial &material,
                                  const Vec3 &normal,
                                  const Vec3 &wo,
                                  unsigned int &rngState,
                                  BsdfSample *outSample)
{
    if (!outSample)
    {
        return false;
    }

    const Vec3 f0 = ComputeF0(material);
    const float specProbRaw = ComputeSpecularProbability(material, f0);
    const float specProb =
        specProbRaw <= 0.0f
            ? 0.0f
            : (specProbRaw < 0.05f ? 0.05f : (specProbRaw > 0.95f ? 0.95f : specProbRaw));

    Vec3 tangent = {};
    Vec3 bitangent = {};
    BuildOrthonormalBasis(normal, tangent, bitangent);

    Vec3 wi = {};
    if (Random01(rngState) < specProb)
    {
        const Vec3 localHalf =
            SampleGgxHalfVector(ComputeAlpha(material), Random01(rngState), Random01(rngState));
        const Vec3 halfVec = Normalize(LocalToWorld(localHalf, tangent, bitangent, normal));
        wi = Normalize(Reflect(Mul(wo, -1.0f), halfVec));
    }
    else
    {
        const Vec3 local = SampleCosineHemisphere(Random01(rngState), Random01(rngState));
        wi = Normalize(LocalToWorld(local, tangent, bitangent, normal));
    }

    if (Dot(normal, wi) <= 0.0f)
    {
        return false;
    }

    float pdf = 0.0f;
    const Vec3 f = EvaluateBsdf(material, normal, wo, wi, &pdf);
    if (pdf <= 1.0e-8f)
    {
        return false;
    }

    outSample->wi = wi;
    outSample->f = f;
    outSample->pdf = pdf;
    outSample->skipNeeMis = ShouldSkipNeeMis(material);
    return true;
}

} // namespace integrator
} // namespace render
} // namespace ybi
