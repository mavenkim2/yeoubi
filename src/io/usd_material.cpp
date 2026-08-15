#include "io/usd_material.h"
#include "scene/openpbr_parameters.h"
#include <pxr/base/gf/vec3f.h>
#include <pxr/base/vt/value.h>
#include <pxr/usd/sdr/shaderNode.h>
#include <pxr/usd/sdr/shaderProperty.h>
#include <pxr/usd/usdShade/material.h>
#include <pxr/usd/usdShade/output.h>

namespace ybi
{

static const pxr::TfToken mtlxToken("mtlx", pxr::TfToken::Immortal);
static const pxr::TfToken mtlxOpenPbr("ND_open_pbr_surface_surfaceshader", pxr::TfToken::Immortal);

struct OpenPbrParameterDescriptor
{
    pxr::TfToken token;

    bool (*assignConstant)(OpenPbrMaterial &, const pxr::VtValue &);

    void (*assignTexture)(OpenPbrMaterial &, uint32_t);
};

static bool AssignConstant(OpenPbrInput<float> &input, const pxr::VtValue &value)
{
    if (!value.IsHolding<float>())
        return false;

    input.constant = value.UncheckedGet<float>();
    input.textureIndex = kInvalidOpenPbrTextureIndex;
    input.source = OpenPbrInputSource::Constant;
    return true;
}

static bool AssignConstant(OpenPbrInput<bool> &input, const pxr::VtValue &value)
{
    if (!value.IsHolding<bool>())
        return false;

    input.constant = value.UncheckedGet<bool>();
    input.textureIndex = kInvalidOpenPbrTextureIndex;
    input.source = OpenPbrInputSource::Constant;
    return true;
}

static bool AssignConstant(OpenPbrInput<Vec3f> &input, const pxr::VtValue &value)
{
    if (!value.IsHolding<pxr::GfVec3f>())
        return false;

    const pxr::GfVec3f &v = value.UncheckedGet<pxr::GfVec3f>();
    input.constant = Vec3f(v[0], v[1], v[2]);
    input.textureIndex = kInvalidOpenPbrTextureIndex;
    input.source = OpenPbrInputSource::Constant;
    return true;
}

template <typename T>
static void AssignTexture(OpenPbrInput<T> &input, uint32_t textureIndex)
{
    input.textureIndex = textureIndex;
    input.source = OpenPbrInputSource::Texture;
}

static const OpenPbrParameterDescriptor openPbrParameterDescriptors[] = {
#define MAKE_PARAMETER_DESCRIPTOR(usdName, type) \
    { \
        pxr::TfToken(#usdName, pxr::TfToken::Immortal), \
        [](OpenPbrMaterial &material, const pxr::VtValue &value) { \
            return AssignConstant(material.usdName, value); \
        }, \
        [](OpenPbrMaterial &material, uint32_t textureIndex) { \
            AssignTexture(material.usdName, textureIndex); \
        }, \
    },
    OPENPBR_PARAMETERS(MAKE_PARAMETER_DESCRIPTOR)
#undef MAKE_PARAMETER_DESCRIPTOR
};

pxr::UsdShadeShader GetMaterialSurfaceShader(const pxr::UsdShadeMaterial &material)
{
    pxr::UsdShadeShader shader;
    pxr::UsdShadeOutput shadeOutput = material.GetSurfaceOutput(mtlxToken);
    if (shadeOutput)
    {
        pxr::UsdShadeAttributeVector valueAttrs =
            pxr::UsdShadeUtils::GetValueProducingAttributes(shadeOutput, true);
        shader = pxr::UsdShadeShader(valueAttrs[0].GetPrim());
    }
    else
    {
        shader = material.ComputeSurfaceSource(pxr::UsdShadeTokens->universalRenderContext);
    }

    return shader;
}

static void AssignDefinitionFallback(const pxr::SdrShaderNodeConstPtr &definition,
                                     const OpenPbrParameterDescriptor &descriptor,
                                     OpenPbrMaterial *material)
{
    if (!definition)
        return;

    const pxr::SdrShaderPropertyConstPtr property = definition->GetShaderInput(descriptor.token);
    if (property)
        descriptor.assignConstant(*material, property->GetDefaultValueAsSdfType());
}

bool LoadOpenPbrMaterial(const pxr::UsdShadeShader &shader,
                         pxr::UsdTimeCode timeCode,
                         OpenPbrMaterial *material)
{
    if (!shader || !material)
        return false;

    pxr::TfToken idToken;
    if (!shader.GetShaderId(&idToken) || idToken != mtlxOpenPbr)
        return false;

    *material = OpenPbrMaterial{};
    const pxr::SdrShaderNodeConstPtr definition = shader.GetShaderNodeForSourceType(mtlxToken);

    for (const OpenPbrParameterDescriptor &descriptor : openPbrParameterDescriptors)
    {
        AssignDefinitionFallback(definition, descriptor, material);

        const pxr::UsdShadeInput input = shader.GetInput(descriptor.token);
        if (!input)
            continue;

        const pxr::UsdShadeAttributeVector producers = input.GetValueProducingAttributes();
        if (producers.empty())
        {
            pxr::VtValue value;
            if (input.Get(&value, timeCode))
                descriptor.assignConstant(*material, value);

            continue;
        }

        for (const pxr::UsdAttribute &producer : producers)
        {
            pxr::VtValue value;
            if (producer.GetPrim().IsA<pxr::UsdShadeShader>())
            {
                pxr::UsdShadeShader sourceShader(producer.GetPrim());
                pxr::TfToken shaderId;
                sourceShader.GetShaderId(&shaderId);
                printf("%s\n", shaderId.GetText());
            }
            else if (producer.Get(&value, timeCode) && descriptor.assignConstant(*material, value))
            {
                break;
            }
        }
    }

    return true;
}

bool SetOpenPbrTexture(OpenPbrMaterial *material,
                       const pxr::TfToken &parameterName,
                       uint32_t textureIndex)
{
    if (!material || textureIndex == kInvalidOpenPbrTextureIndex)
        return false;

    for (const OpenPbrParameterDescriptor &descriptor : openPbrParameterDescriptors)
    {
        if (descriptor.token == parameterName)
        {
            descriptor.assignTexture(*material, textureIndex);
            return true;
        }
    }

    return false;
}

} // namespace ybi
