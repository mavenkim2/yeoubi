#pragma once

#include "render/integrator_common.h"
#include "render/integrator_shading.h"
#include "render/openpbr_ybi.h"
#include "scene/material_light.h"
#include "util/math_common.h"

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
    Vec3 specularColor = Vec3(1.0f, 1.0f, 1.0f);
    float roughness = 1.0f;
    float metallic = 0.0f;
    float ior = 1.5f;
    float opacity = 1.0f;
    float clearcoat = 0.0f;
    float clearcoatRoughness = 0.01f;
    float opacityThreshold = 0.0f;
    unsigned int flags = 0u;
    unsigned int useSpecularWorkflow = 0u;
};

struct PreparedMaterial
{
    OpenPBR_PreparedBsdf prepared = {};
    Vec3 geomNormal = Vec3(0.0f, 0.0f, 1.0f);
    Vec3 shadingNormal = Vec3(0.0f, 0.0f, 1.0f);
    Vec3 tangent = Vec3(1.0f, 0.0f, 0.0f);
    Vec3 bitangent = Vec3(0.0f, 1.0f, 0.0f);
    Vec3 emissiveColor = Vec3(0.0f, 0.0f, 0.0f);
    float opacity = 1.0f;
    float opacityThreshold = 0.0f;
    unsigned int flags = 0u;
};

struct BsdfSample
{
    Vec3 wi = Vec3(0.0f, 0.0f, 1.0f);
    Vec3 weight = Vec3(0.0f, 0.0f, 0.0f);
    float pdf = 0.0f;
    bool skipNeeMis = false;
    bool transmission = false;
};

YBI_INTEGRATOR_HD float MaxFloat(float a, float b)
{
    return a > b ? a : b;
}

YBI_INTEGRATOR_HD float MinFloat(float a, float b)
{
    return a < b ? a : b;
}

YBI_INTEGRATOR_HD float SafeSqrt(float x)
{
    return sqrtf(MaxFloat(x, 0.0f));
}

YBI_INTEGRATOR_HD Vec3 ClampVec3(const Vec3 &v, float lo, float hi)
{
    return Clamp(v, lo, hi);
}

YBI_INTEGRATOR_HD bool ShouldAlphaCutout(const PreparedMaterial &material)
{
    return material.opacityThreshold > 0.0f && material.opacity < material.opacityThreshold;
}

YBI_INTEGRATOR_HD OpenPBR_Basis MakeOpenPbrBasis(const ShadingFrame &frame)
{
    OpenPBR_Basis basis = {};
    basis.t = frame.tangent;
    basis.b = frame.bitangent;
    basis.n = frame.shadingNormal;
    return basis;
}

YBI_INTEGRATOR_HD void FillOpenPbrEmission(const Vec3 &emissiveColor,
                                           OpenPBR_ResolvedInputs *inputs)
{
    if (!inputs)
    {
        return;
    }

    const float luminance = MaxFloat(Luminance(emissiveColor), 0.0f);
    inputs->emission_luminance = luminance;
    inputs->emission_color =
        luminance > 1.0e-6f ? emissiveColor / luminance : Vec3(1.0f, 1.0f, 1.0f);
}

YBI_INTEGRATOR_HD OpenPBR_ResolvedInputs BuildOpenPbrResolvedInputs(const EvaluatedMaterial &material,
                                                                    const ShadingFrame &frame)
{
    OpenPBR_ResolvedInputs inputs = openpbr_make_default_resolved_inputs();
    inputs.base_color = material.baseColor;
    inputs.base_metalness = material.useSpecularWorkflow != 0u ? 0.0f : material.metallic;
    inputs.specular_color =
        material.useSpecularWorkflow != 0u ? material.specularColor : Vec3(1.0f, 1.0f, 1.0f);
    inputs.specular_roughness = material.roughness;
    inputs.specular_ior = material.ior;
    inputs.coat_weight = material.clearcoat;
    inputs.coat_roughness = material.clearcoatRoughness;
    inputs.coat_color = Vec3(1.0f, 1.0f, 1.0f);
    inputs.coat_ior = 1.5f;
    inputs.transmission_weight = material.opacity < 1.0f ? (1.0f - material.opacity) : 0.0f;
    inputs.transmission_color = Vec3(1.0f, 1.0f, 1.0f);
    inputs.transmission_depth = 0.0f;
    inputs.geometry_opacity = material.opacity;
    inputs.geometry_thin_walled = inputs.transmission_weight > 0.0f;
    inputs.shading_basis = MakeOpenPbrBasis(frame);
    inputs.coat_basis = inputs.shading_basis;
    FillOpenPbrEmission(material.emissiveColor, &inputs);
    return inputs;
}

YBI_INTEGRATOR_HD PreparedMaterial PrepareMaterialBsdf(const EvaluatedMaterial &material,
                                                       const ShadingFrame &frame,
                                                       const Vec3 &pathThroughput,
                                                       const Vec3 &viewDirection)
{
    PreparedMaterial out = {};
    out.geomNormal = frame.geomNormal;
    out.shadingNormal = frame.shadingNormal;
    out.tangent = frame.tangent;
    out.bitangent = frame.bitangent;
    out.emissiveColor = material.emissiveColor;
    out.opacity = material.opacity;
    out.opacityThreshold = material.opacityThreshold;
    out.flags = material.flags;

    const OpenPBR_ResolvedInputs inputs = BuildOpenPbrResolvedInputs(material, frame);
    out.prepared = openpbr_prepare_bsdf_and_volume(
        inputs, pathThroughput, OpenPbrDefaultRgbWavelengthsNm(), 1.0f, viewDirection);
    return out;
}

YBI_INTEGRATOR_HD Vec3 EvaluateBsdf(const PreparedMaterial &material,
                                    const Vec3 &wi,
                                    float *outPdf = nullptr)
{
    if (outPdf)
    {
        *outPdf = openpbr_pdf(material.prepared, wi);
    }
    return OpenPbrCombinedWeight(openpbr_eval(material.prepared, wi));
}

YBI_INTEGRATOR_HD bool SampleBsdf(const PreparedMaterial &material,
                                  unsigned int &rngState,
                                  BsdfSample *outSample)
{
    if (!outSample)
    {
        return false;
    }

    vec3 lightDirection = vec3(0.0f, 0.0f, 1.0f);
    OpenPBR_DiffuseSpecular weight = openpbr_make_zero_diffuse_specular();
    float pdf = 0.0f;
    OpenPBR_BsdfLobeType sampledType = OpenPBR_BsdfLobeTypeNone;
    openpbr_sample(material.prepared,
                   vec3(Random01(rngState), Random01(rngState), Random01(rngState)),
                   lightDirection,
                   weight,
                   pdf,
                   sampledType);
    if (pdf <= 1.0e-8f)
    {
        return false;
    }

    outSample->wi = lightDirection;
    outSample->weight = OpenPbrCombinedWeight(weight);
    outSample->pdf = pdf;
    outSample->skipNeeMis = OpenPbrSampleIsDelta(sampledType);
    outSample->transmission = (sampledType & OpenPBR_BsdfLobeTypeTransmission) != 0u;
    return MaxComponent(outSample->weight) > 0.0f;
}

} // namespace integrator
} // namespace render
} // namespace ybi
