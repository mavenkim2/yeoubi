#pragma once

#include "coat.h"
#include "core/math.h"
#include "ggx.h"
#include "kernel/bsdf/openpbr/bsdf_types.h"
#include "kernel/bsdf/openpbr/coat.h"
#include "kernel/bsdf/openpbr/sample_table.h"
#include "kernel/launch/shading_context.h"
#include "kernel/shading/material_shaders.h"
#include "kernel/shading/shading.h"
#include "kernel/scene.h"

namespace ybi
{

YBI_GPU Vec2f ComputeAnisotropicRoughness(float roughness, float anisotropy)
{
    Vec2f anisotropicRoughness;
    anisotropicRoughness.x = Sqr(roughness) * sqrtf(2.f / (1.f + Sqr(1.0f - anisotropy)));
    anisotropicRoughness.y = (1.0f - anisotropy) * anisotropicRoughness.x;

    // TODO: support dirac delta lobes
    anisotropicRoughness = Max(anisotropicRoughness, Sqr(0.001f));
    return anisotropicRoughness;
}

YBI_GPU Vec3f CalculateThinWallTransmissionColor(const Vec3f &viewDirection,
                                                 float specularIor,
                                                 const Vec3f &transmissionColor)
{
    float cosTheta = Abs(viewDirection.z);
    float sin2Theta_i = 1.f - Sqr(cosTheta);

    // TODO: this ior has to change when the coat layer is added
    float sin2Theta_t = sin2Theta_i / Sqr(specularIor);

    if (sin2Theta_t >= 1.f)
    {
        return Vec3f(0.f);
    }

    float cos_t = sqrtf(1.f - sin2Theta_t);
    float pathLength = 1.f / cos_t;
    Vec3f thinWallTransmissionColor = Pow(transmissionColor, pathLength);
    return thinWallTransmissionColor;
}

// https://blog.selfshadow.com/publications/s2017-shading-course/imageworks/s2017_pbs_imageworks_slides_v2.pdf,
// slide 18.
YBI_GPU float ApproximateAverageFresnel(const float eta)
{
    if (eta > 1.0f)
        return (eta - 1.0f) / (4.08567f + 1.00071f * eta);
    else
        return 0.997118f + 0.1014f * eta - 0.965241f * (eta * eta) - 0.130607f * (eta * eta) * eta;
}

// NOTE: ior weighting is only used for reflection, it is not applied to refraction
YBI_GPU float ApplyWeightToIor(float ior, float weight)
{
    float fNormalIncidence = Sqr((1 - ior) / (1 + ior));
    float eps = sqrtf(Min(weight * fNormalIncidence, 0.9999f));
    eps = ior >= 1.f ? eps : -eps;
    float result = (1 + eps) / (1 - eps);
    return result;
}

YBI_GPU bool InitializeOpenPbrBsdf(const ShadingContext &context,
                                   const DeviceMaterial2 &material,
                                   const ShadingInfo &info,
                                   OpenPbrBsdf &bsdf,
                                   const Vec3f &worldNormal,
                                   const Vec3f &worldViewDirection,
                                   OrthonormalBasis tangentFrame,
                                   const Vec3f &throughput,
                                   float u)
{
    const float geometryOpacity = EvaluateShaderFloat(context, info, material.geometry_opacity);
    if (u >= geometryOpacity)
    {
        return false;
    }

    const float fuzzWeight = EvaluateShaderFloat(context, info, material.fuzz_weight);
    const Vec3f fuzzColor = EvaluateShaderVec3(context, info, material.fuzz_color);
    const float fuzzRoughness = EvaluateShaderFloat(context, info, material.fuzz_roughness);

#if 0
    const float geometryOpacity = 1.f;
    if (u >= geometryOpacity)
    {
        return false;
    }

    const float fuzzWeight = 0.f;
    const Vec3f fuzzColor = Vec3f(0.f);
    const float fuzzRoughness = 0.f;

    const float specularIor = 1.5f;
    const float specularWeight = 1.f;
    const float specularAnisotropy = 0.f;
    const float specularRoughness = 0.3f;
    const Vec3f specularColor = Vec3f(1.f);

    const Vec3f baseColor = Vec3f(0.8f);
    const float baseWeight = 1.f;

    // Translucent lobe
    const Vec3f transmissionColor = Vec3f(.5f, 1.f, .5f);
    const float transmissionWeight = .9f;

    // Coat
    const Vec3f coatColor = Vec3f(1.f);
    const float coatWeight = 0.f;
    const float coatRoughness = 0.f;
    const float coatRoughnessAnisotropy = 0.f;
    const float coatDarkening = 1.f;
    const float coatIor = 1.6f;

    // Thin wall
    const bool isThinWall = true;
#endif

    // const float fuzzWeight = 0.f;
    // const Vec3f fuzzColor = Vec3f(1.f);
    // const float fuzzRoughness = 0.f;
    //
    // const float specularIor = 1.5f;
    // const float specularWeight = 0.f; // 1.f;
    // const float specularAnisotropy = 0.f;
    // const float specularRoughness = 0.3f;
    // const Vec3f specularColor = Vec3f(1.f);

    const float specularIor = EvaluateShaderFloat(context, info, material.specular_ior);
    const float specularWeight = EvaluateShaderFloat(context, info, material.specular_weight);
    const float specularAnisotropy =
        EvaluateShaderFloat(context, info, material.specular_roughness_anisotropy);
    const float specularRoughness = EvaluateShaderFloat(context, info, material.specular_roughness);
    const Vec3f specularColor = EvaluateShaderVec3(context, info, material.specular_color);

    // const Vec3f baseColor = Vec3f(0.8f);
    // const float baseWeight = 1.f;
    const Vec3f baseColor = EvaluateShaderVec3(context, info, material.base_color);
    const float baseWeight = EvaluateShaderFloat(context, info, material.base_weight);

    // Translucent lobe
    const Vec3f transmissionColor = EvaluateShaderVec3(context, info, material.transmission_color);
    const float transmissionWeight =
        EvaluateShaderFloat(context, info, material.transmission_weight);

    // Coat
    const Vec3f coatColor = EvaluateShaderVec3(context, info, material.coat_color);
    const float coatWeight = EvaluateShaderFloat(context, info, material.coat_weight);
    const float coatRoughness = EvaluateShaderFloat(context, info, material.coat_roughness);
    const float coatRoughnessAnisotropy =
        EvaluateShaderFloat(context, info, material.coat_roughness_anisotropy);
    const float coatDarkening = EvaluateShaderFloat(context, info, material.coat_darkening);
    const float coatIor = EvaluateShaderFloat(context, info, material.coat_ior);
    // const Vec3f coatColor = Vec3f(1.f);
    // const float coatWeight = 0.f;
    // const float coatRoughness = 0.f;
    // const float coatRoughnessAnisotropy = 0.f;
    // const float coatDarkening = 1.f;
    // const float coatIor = 1.6f;

    // Thin wall
    const bool isThinWall = EvaluateShaderBool(context, info, material.geometry_thin_walled);

    bsdf.fuzzWeight = fuzzWeight;
    bsdf.fuzzColor = fuzzColor;
    bsdf.fuzzRoughness = fuzzRoughness;

    const bool backFacing = Dot(worldViewDirection, worldNormal) < 0.f;
    if (backFacing)
    {
        tangentFrame.x = -tangentFrame.x;
        tangentFrame.y = -tangentFrame.y;
        tangentFrame.z = -tangentFrame.z;
    }
    bsdf.tangentFrame = tangentFrame;
    const Vec3f localViewDirection = tangentFrame.WorldToLocal(worldViewDirection);

    bsdf.specularColor = specularColor;
    bsdf.baseColor = baseColor * baseWeight;

    float opaqueDielectricWeight = 1.f - transmissionWeight;
    bsdf.baseColor *= opaqueDielectricWeight;

    bool inside = backFacing && !isThinWall;
    bsdf.inside = inside;

    // SampleOpaqueDielectricTable(params.opaqueDielectricEnergyComplementTexture,
    //                             float ior,
    //                             viewDirection.z,
    //                             specularRoughness);

    // Coat roughening
    // TODO: handle nested dielectrics
    float exteriorIor = 1.f;
    float relativeCoatIor = 1.f;
    bool enterCoat = coatWeight > 0.f && !inside;
    if (enterCoat)
    {
        // TODO: why doesn't this used the weighted coat ior
        relativeCoatIor = coatIor / exteriorIor;
    }

    // Calculate iors
    float iorForTransReflection;
    float iorForRefraction;
    float iorForOpaqueReflection;

    {
        // const float relativeIorForTransReflection = Lerp();
        const float weightedCoatIor =
            enterCoat ? Lerp(exteriorIor, coatIor, coatWeight) : exteriorIor;
        float weightedSpecularIor =
            weightedCoatIor * ApplyWeightToIor(specularIor / weightedCoatIor, specularWeight);

        if (enterCoat)
        {
            iorForTransReflection = weightedSpecularIor / weightedCoatIor;
            if (iorForTransReflection < 1.f)
            {
                iorForTransReflection = 1.f / iorForTransReflection;
            }
            iorForRefraction = weightedSpecularIor / exteriorIor;
            iorForOpaqueReflection = iorForTransReflection;
        }
        else
        {
            iorForTransReflection = exteriorIor / weightedSpecularIor;
            iorForRefraction = iorForTransReflection;
            iorForOpaqueReflection = weightedSpecularIor / exteriorIor;
        }

        const float MinRefractionIorAboveOne = 1.0001f;
        if (iorForRefraction >= 1.0f)
        {
            iorForRefraction = Max(iorForRefraction, MinRefractionIorAboveOne);
        }
        else
        {
            const float InverseMinRefractionIorAboveOne = 1.0f / MinRefractionIorAboveOne;
            iorForRefraction = Min(iorForRefraction, InverseMinRefractionIorAboveOne);
        }

        bsdf.specularReflectionIor = iorForOpaqueReflection;
        bsdf.specularRefractionIor = iorForRefraction;
    }

    // Heuristic roughness adjustments
    float adjustedSpecularRoughness;
    float adjustedCoatRoughness;
    {
        // Adjust coat roughness to account for fuzz roughness
        const float fuzzEffectOnUnderlyingRoughness = 0.005f;
        float avgFuzzColor = (fuzzColor.x + fuzzColor.y + fuzzColor.z) / 3.f;
        const float fuzzRoughnessAdjustmentFactor =
            avgFuzzColor * fuzzRoughness * fuzzEffectOnUnderlyingRoughness;
        adjustedCoatRoughness = Lerp(
            coatRoughness,
            Root4(Min(1.f,
                      Pow4(coatRoughness) + fuzzRoughnessAdjustmentFactor * Pow4(fuzzRoughness))),
            fuzzWeight);

        // Adjust specular roughness to account for coat roughness
        const float coatRoughnessAdjustmentFactor =
            Lerp(1.f, 0.f, relativeCoatIor >= 1.f ? 1.f / relativeCoatIor : relativeCoatIor);

        adjustedSpecularRoughness =
            Lerp(specularRoughness,
                 Root4(Min(1.f,
                           Pow4(specularRoughness) +
                               coatRoughnessAdjustmentFactor * Pow4(adjustedCoatRoughness))),
                 coatWeight);

        // Calculate roughness
        bsdf.specularRoughness =
            ComputeAnisotropicRoughness(adjustedSpecularRoughness, specularAnisotropy);
        bsdf.coatRoughness =
            ComputeAnisotropicRoughness(adjustedCoatRoughness, coatRoughnessAnisotropy);
    }

    // Calculate coat darkening
    // See section 3.4.1 in OpenPBR spec
    {
        // TODO: SSS, metal
        // NOTE: openpbr approximates this with the average fresnel
        float Ks = FrDielectric(localViewDirection.z, relativeCoatIor);
        float Kr = 1.f - (1.f - ApproximateAverageFresnel(relativeCoatIor)) / Sqr(relativeCoatIor);

        const float specularityOpaqueDielectric =
            FrDielectric(localViewDirection.z, iorForOpaqueReflection);
        const float specularitySSSTrans =
            FrDielectric(localViewDirection.z, iorForTransReflection);
        // const float specularity_dielectric = opaque_dielectric * specularity_opaque_dielectric +
        //                                      (sss + trans) * specularity_sss_trans;
        const float specularityDielectric = opaqueDielectricWeight * specularityOpaqueDielectric +
                                            transmissionWeight * specularitySSSTrans;

        const float specularityMetal = 0.f; // TODO
        const float specularityBase = specularityDielectric + specularityMetal;
        const float effectiveRoughnessBase = Lerp(1.f, adjustedSpecularRoughness, specularityBase);

        float K = Lerp(Ks, Kr, effectiveRoughnessBase);

        Vec3f baseAlbedoFromOpaqueDielectric =
            opaqueDielectricWeight *
            Lerp(baseColor * baseWeight, 1.f, specularityOpaqueDielectric);
        Vec3f baseAlbedoFromTrans = Vec3f(transmissionWeight);

        Vec3f Eb = baseAlbedoFromTrans + baseAlbedoFromOpaqueDielectric;
        Vec3f delta = (1.f - K) / (1.f - Eb * K);
        Vec3f uniformModulatedCoatDarkeningFactor = Lerp(1.f, delta, coatWeight * coatDarkening);

        const float coatEnergyReflected =
            CalculateCoatEnergyReflected(context.opaqueDielectricEnergyComplementTexture,
                                         relativeCoatIor,
                                         localViewDirection.z,
                                         adjustedCoatRoughness,
                                         coatWeight);
        bsdf.coatColor = coatColor;
        bsdf.coatWeight = coatWeight;
        bsdf.coatEnergyReflected = coatEnergyReflected;
        bsdf.coatDarkeningFactor =
            uniformModulatedCoatDarkeningFactor *
            ScaleCoatColor(coatColor, localViewDirection.z, relativeCoatIor, coatWeight) *
            (1.f - coatEnergyReflected);
        bsdf.coatIor = relativeCoatIor;
    }

    // Thin wall
    {
        bsdf.thinWall = isThinWall;
        if (isThinWall)
        {
            // NOTE: the transmission color represents the transmittance through the thin
            // sheet at normal incidence.
            Vec3f thinWallTransmissionColor = CalculateThinWallTransmissionColor(
                localViewDirection, bsdf.specularRefractionIor, transmissionColor);

            // Calculate the total transmittance through the thin sheet.
            // T = (1-R)^2 * (1 + R^2 + R^4 ...) = (1 - R) / (1 + R)
            float thinWallFresnel =
                FrDielectric(Abs(localViewDirection.z), bsdf.specularRefractionIor);
            float totalThinWallTransmittance = (1.f - thinWallFresnel) / (1.f + thinWallFresnel);
            thinWallTransmissionColor *= transmissionWeight * totalThinWallTransmittance;

            bsdf.thinWallTransmissionColor = thinWallTransmissionColor;
        }
        else
        {
            bsdf.thinWallTransmissionColor = Vec3f(0.f);
        }
    }

    // Calculate lobe weights
    // Thin Wall Transmission Lobe
    bsdf.lobeWeights[THIN_WALL_TRANS] =
        (throughput * bsdf.thinWallTransmissionColor).MaxComponent();
    bsdf.lobeWeights[DIFFUSE] = opaqueDielectricWeight;
    bsdf.lobeWeights[SPECULAR] = 0.5f;
    return true;
}

} // namespace ybi
