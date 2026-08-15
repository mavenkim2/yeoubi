#include "io/usd_material.h"
#include <pxr/base/gf/vec3f.h>
#include <pxr/usd/sdr/shaderNode.h>
#include <pxr/usd/sdr/shaderProperty.h>
#include <pxr/usd/usdShade/material.h>
#include <pxr/usd/usdShade/nodeGraph.h>
#include <pxr/usd/usdShade/output.h>

namespace ybi
{

static const pxr::TfToken mtlxToken("mtlx", pxr::TfToken::Immortal);
static const pxr::TfToken previewSurfaceToken("UsdPreviewSurface", pxr::TfToken::Immortal);
static const pxr::TfToken mtlxOpenPbr("ND_open_pbr_surface_surfaceshader", pxr::TfToken::Immortal);

static const pxr::TfToken openPbrParameterTokens[] = {
    pxr::TfToken("base_weight", pxr::TfToken::Immortal),
    pxr::TfToken("base_color", pxr::TfToken::Immortal),
    pxr::TfToken("base_metalness", pxr::TfToken::Immortal),
    pxr::TfToken("base_diffuse_roughness", pxr::TfToken::Immortal),

    pxr::TfToken("specular_weight", pxr::TfToken::Immortal),
    pxr::TfToken("specular_color", pxr::TfToken::Immortal),
    pxr::TfToken("specular_roughness", pxr::TfToken::Immortal),
    pxr::TfToken("specular_roughness_anisotropy", pxr::TfToken::Immortal),
    pxr::TfToken("specular_ior", pxr::TfToken::Immortal),

    pxr::TfToken("transmission_weight", pxr::TfToken::Immortal),
    pxr::TfToken("transmission_color", pxr::TfToken::Immortal),
    pxr::TfToken("transmission_depth", pxr::TfToken::Immortal),
    pxr::TfToken("transmission_scatter", pxr::TfToken::Immortal),
    pxr::TfToken("transmission_scatter_anisotropy", pxr::TfToken::Immortal),
    pxr::TfToken("transmission_dispersion_scale", pxr::TfToken::Immortal),
    pxr::TfToken("transmission_dispersion_abbe_number", pxr::TfToken::Immortal),

    pxr::TfToken("subsurface_weight", pxr::TfToken::Immortal),
    pxr::TfToken("subsurface_color", pxr::TfToken::Immortal),
    pxr::TfToken("subsurface_radius", pxr::TfToken::Immortal),
    pxr::TfToken("subsurface_radius_scale", pxr::TfToken::Immortal),
    pxr::TfToken("subsurface_scatter_anisotropy", pxr::TfToken::Immortal),

    pxr::TfToken("coat_weight", pxr::TfToken::Immortal),
    pxr::TfToken("coat_color", pxr::TfToken::Immortal),
    pxr::TfToken("coat_roughness", pxr::TfToken::Immortal),
    pxr::TfToken("coat_roughness_anisotropy", pxr::TfToken::Immortal),
    pxr::TfToken("coat_ior", pxr::TfToken::Immortal),
    pxr::TfToken("coat_darkening", pxr::TfToken::Immortal),

    pxr::TfToken("fuzz_weight", pxr::TfToken::Immortal),
    pxr::TfToken("fuzz_color", pxr::TfToken::Immortal),
    pxr::TfToken("fuzz_roughness", pxr::TfToken::Immortal),

    pxr::TfToken("emission_luminance", pxr::TfToken::Immortal),
    pxr::TfToken("emission_color", pxr::TfToken::Immortal),

    pxr::TfToken("thin_film_weight", pxr::TfToken::Immortal),
    pxr::TfToken("thin_film_thickness", pxr::TfToken::Immortal),
    pxr::TfToken("thin_film_ior", pxr::TfToken::Immortal),

    pxr::TfToken("geometry_opacity", pxr::TfToken::Immortal),
    pxr::TfToken("geometry_thin_walled", pxr::TfToken::Immortal),
    pxr::TfToken("geometry_normal", pxr::TfToken::Immortal),
    pxr::TfToken("geometry_tangent", pxr::TfToken::Immortal),
    pxr::TfToken("geometry_coat_normal", pxr::TfToken::Immortal),
    pxr::TfToken("geometry_coat_tangent", pxr::TfToken::Immortal),
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

static void ProcessOpenPBRMaterial(const pxr::UsdShadeShader &shader, pxr::UsdTimeCode timeCode)
{
    const size_t numParameters =
        sizeof(openPbrParameterTokens) / sizeof(openPbrParameterTokens[0]);

    for (size_t i = 0; i < numParameters; i++)
    {
        if (pxr::UsdShadeInput input = shader.GetInput(openPbrParameterTokens[i]))
        {
            auto sources = input.GetConnectedSources();

            // for (auto &source : sources)
            // {
            //     pxr::UsdPrim sourcePrim = source.source.GetPrim();
            //     if (sourcePrim.IsA<pxr::UsdShadeShader>())
            //     {
            //         pxr::UsdShadeShader sourceShader(sourcePrim);
            //
            //         pxr::TfToken shaderId;
            //         sourceShader.GetShaderId(&shaderId);
            //
            //         printf("%s\n", shaderId.GetText());
            //     }
            //     else if (sourcePrim.IsA<pxr::UsdShadeNodeGraph>())
            //     {
            //         pxr::UsdShadeNodeGraph nodeGraph(sourcePrim);
            //         sourcePrim.pxr::UsdShadeShader sourceShader(sourcePrim);
            //
            //         pxr::TfToken shaderId;
            //         sourceShader.GetShaderId(&shaderId);
            //         printf("sigh\n");
            //     }
            //     else
            //     {
            //         printf("why\n");
            //     }
            // }

            const pxr::UsdShadeAttributeVector producers = input.GetValueProducingAttributes();
            if (producers.size())
            {
                for (const pxr::UsdAttribute &producer : producers)
                {
                    printf("%s\n", producer.GetName().GetText());
                    pxr::VtValue value;
                    if (producer.GetPrim().IsA<pxr::UsdShadeShader>())
                    {
                        pxr::UsdShadeShader sourceShader(producer.GetPrim());
                        pxr::TfToken shaderId;
                        sourceShader.GetShaderId(&shaderId);
                        printf("%s\n", shaderId.GetText());
                    }
                    else if (producer.Get(&value, timeCode))
                    {
                        if (value.IsHolding<float>())
                        {
                            if (i == 21)
                            {
                                printf("%f\n", value.UncheckedGet<float>());
                            }
                        }
                        else if (value.IsHolding<pxr::GfVec3f>())
                        {
                            pxr::GfVec3f v = value.UncheckedGet<pxr::GfVec3f>();
                            // printf("%f %f %f\n", v[0], v[1], v[2]);
                        }
                    }
                }
            }
            else
            {
                pxr::SdrShaderNodeConstPtr definition =
                    shader.GetShaderNodeForSourceType(mtlxToken);

                if (definition)
                {
                    pxr::SdrShaderPropertyConstPtr property =
                        definition->GetShaderInput(input.GetBaseName());

                    if (property)
                    {
                        const pxr::VtValue &defaultValue = property->GetDefaultValueAsSdfType();

                        if (defaultValue.IsHolding<float>())
                        {
                            // printf("%f\n", defaultValue.UncheckedGet<float>());
                        }
                        else if (defaultValue.IsHolding<pxr::GfVec3f>())
                        {
                            pxr::GfVec3f v = defaultValue.UncheckedGet<pxr::GfVec3f>();
                            // printf("%f %f %f\n", v[0], v[1], v[2]);
                        }
                    }
                }
            }
        }
    }
    printf("\n");
}

void ProcessUsdShader(const pxr::UsdShadeShader &shader, pxr::UsdTimeCode timeCode)
{
    pxr::TfToken idToken;
    if (shader && shader.GetShaderId(&idToken))
    {
        if (idToken == previewSurfaceToken)
        {
        }
        else if (idToken == mtlxOpenPbr)
        {
            ProcessOpenPBRMaterial(shader, timeCode);
        }
    }
}

} // namespace ybi
