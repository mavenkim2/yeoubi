#pragma once

#if 1
#include "device/cuda_compat.h"

namespace ybi
{

YBI_GPU Vec4f CalculateShadowTransmittance(const KernelParams &params,
                                           const DeviceMaterial2 &material,
                                           const ShadingInfo &info,
                                           const OrthonormalBasis &tangentFrame,
                                           const Vec3f &worldViewDirection,
                                           const Vec3f &worldNormal)
{
    const Vec3f localViewDirection = tangentFrame.WorldToLocal(worldViewDirection);
    const float cosTheta = Dot(worldNormal, worldViewDirection);

    const bool thinWall = EvaluateShader<bool>(params, info, material.geometry_thin_walled);

    const float transmissionWeight =
        EvaluateShader<float>(params, info, material.transmission_weight);
    const float subsurfaceWeight =
        EvaluateShader<float>(params, info, material.subsurface_weight);
    const float metalness = EvaluateShader<float>(params, info, material.base_metalness);
    const float subsurfaceTranslucency = thinWall ? 0.f : subsurfaceWeight;
    const float translucency = Lerp(subsurfaceTranslucency, 1.f, transmissionWeight);

    const float specularRoughness =
        EvaluateShader<float>(params, info, material.specular_roughness);
    const float prob = translucency * (1.f - metalness) * (1.f - specularRoughness);

    if (prob == 0.f)
    {
        return Vec4f(0.f);
    }

    // NOTE: Evaluate only parameters necessary to obtain shadow transmittance
    const Vec3f transmissionColor =
        EvaluateShader<Vec3f>(params, info, material.transmission_color);

    const float specularWeight = EvaluateShader<float>(params, info, material.specular_weight);
    const float specularIor = EvaluateShader<float>(params, info, material.specular_ior);

    const Vec3f coatColor = EvaluateShader<Vec3f>(params, info, material.coat_color);
    const float coatWeight = EvaluateShader<float>(params, info, material.coat_weight);

    const float thinWallAwareCosTheta = thinWall ? Abs(cosTheta) : cosTheta;

    const float weightedSpecularIor = ApplyWeightToIor(specularIor, specularWeight);
    const float externalIor =
        weightedSpecularIor >= 1.0f ? weightedSpecularIor : 1.0f / weightedSpecularIor;

    float fresnelFactor;
    if (thinWallAwareCosTheta > 0.0f)
    {
        const float fresnelSmooth = FrDielectric(thinWallAwareCosTheta, externalIor);
        const float fresnelRough = ApproximateAverageFresnel(externalIor);
        const float fresnel = Lerp(fresnelSmooth, fresnelRough, specularRoughness);
        const float fresnelTransmissionOneSurface = 1.0f - fresnel;
        fresnelFactor = Sqr(fresnelTransmissionOneSurface);
    }
    else
    {
        fresnelFactor = 1.0f;
    }

    Vec3f finalTransmission = Vec3f(fresnelFactor);

    if (coatWeight > 0.0f && coatColor.MinComponent() < 1.0f)
    {
        const Vec3f fractionalCoatColor = Lerp(1.f, coatColor, coatWeight);

        const Vec3f sqrtFractionalCoatColor = Sqrt(fractionalCoatColor);
        finalTransmission *= sqrtFractionalCoatColor;
    }

    const float transmissionDepth =
        EvaluateShader<float>(params, info, material.transmission_depth);
    const Vec3f transmissionTint =
        thinWall ? CalculateThinWallTransmissionColor(
                       localViewDirection, externalIor, transmissionColor)
                 : transmissionDepth == 0.f ? transmissionColor : Vec3f(1.f);

    finalTransmission *= Lerp(1.f, transmissionTint, transmissionWeight / translucency);

    return Vec4f(finalTransmission, prob);
}
} // namespace ybi
#endif

#if 0
#include "kernel/bsdf/openpbr/resolved_inputs.h"

namespace ybi
{

YBI_GPU Vec4f CalculateShadowTransmittance(const KernelParams &kernelParams,
                                           const DeviceMaterial2 &material,
                                           const ShadingInfo &info,
                                           const OrthonormalBasis &tangentFrame,
                                           const Vec3f &worldViewDirection,
                                           const Vec3f &worldNormal)
{
    const OpenPBR_ResolvedInputs inputs =
        ResolveOpenPbrInputs(kernelParams, material, info, tangentFrame);
    OpenPBR_VolumeDerivedProps volumeDerivedProps;
    OpenPBR_PreparedBsdf prepared;
    openpbr_prepare_volume(inputs, volumeDerivedProps, prepared, true);
    const vec4 result = openpbr_translucent_shadow_weight_and_prob(
        inputs, volumeDerivedProps, Dot(worldNormal, worldViewDirection));
    return Vec4f(result.x, result.y, result.z, result.w);
}

} // namespace ybi
#endif
