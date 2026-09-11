#pragma once

#include "core/math.h"
#include "device/kernel_params.h"
#include "render/shaders/material_shaders.h"
#include "render/shaders/shading.h"
#include "scene/device_scene.h"

namespace ybi
{

struct OpenPbrBsdf
{
    float fuzzWeight;
    Vec3f fuzzColor;
    float fuzzRoughness;

    float specularIor;
    Vec2f specularRoughness;
    Vec3f specularColor;

    Vec3f baseColor;
};

YBI_GPU Vec2f ComputeAnisotropicRoughness(float roughness, float anisotropy)
{
    Vec2f anisotropicRoughness;
    anisotropicRoughness.x = Sqr(roughness) * sqrtf(2.f / (1.f + Sqr(1.0f - anisotropy)));
    anisotropicRoughness.y = (1.0f - anisotropy) * anisotropicRoughness.x;

    // TODO: support dirac delta lobes
    anisotropicRoughness = Max(anisotropicRoughness, Sqr(0.001f));
    return anisotropicRoughness;
}

YBI_GPU void InitializeOpenPbrBsdf(const KernelParams &params,
                                   const DeviceMaterial2 &material,
                                   const ShadingInfo &info,
                                   OpenPbrBsdf &bsdf)
{
    bsdf.fuzzWeight = EvaluateShader<float>(params, info, material.fuzz_weight);
    bsdf.fuzzColor = EvaluateShader<Vec3f>(params, info, material.fuzz_color);
    bsdf.fuzzRoughness = EvaluateShader<float>(params, info, material.fuzz_roughness);

    float specularIor = EvaluateShader<float>(params, info, material.specular_ior);
    float specularWeight = EvaluateShader<float>(params, info, material.specular_weight);
    float specularAnisotropy =
        EvaluateShader<float>(params, info, material.specular_roughness_anisotropy);
    float specularRoughness = EvaluateShader<float>(params, info, material.specular_roughness);
    Vec3f specularColor = EvaluateShader<Vec3f>(params, info, material.specular_color);

    Vec3f baseColor = EvaluateShader<Vec3f>(params, info, material.base_color);
    float baseWeight = EvaluateShader<float>(params, info, material.base_weight);

    // Calculate base dielectric ior
    float fNormalIncidence = Sqr((1 - specularIor) / (1 + specularIor));
    float eps = sqrtf(Min(specularWeight * fNormalIncidence, 0.9999f));
    eps = specularIor >= 1.f ? eps : -eps;
    bsdf.specularIor = (1 + eps) / (1 - eps);

    // Calculate roughness
    bsdf.specularRoughness = ComputeAnisotropicRoughness(specularRoughness, specularAnisotropy);

    bsdf.specularColor = specularColor;
    bsdf.baseColor = baseColor * baseWeight;
}

} // namespace ybi
