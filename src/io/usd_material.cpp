#include "io/usd_material.h"
#include "core/array.h"
#include "core/assert.h"
#include "core/matrix4x4.h"
#include "core/vec2f.h"
#include "core/vec4f.h"
#include "io/usd_helpers.h"
#include "pxr/usd/usdShade/udimUtils.h"
#include "render/material_opcodes.h"
#include "scene/material.h"
#include "scene/openpbr_parameters.h"
#include "scene/scene.h"
#include <memory>
#include <pxr/base/gf/matrix3d.h>
#include <pxr/base/gf/matrix4d.h>
#include <pxr/base/gf/vec2f.h>
#include <pxr/base/gf/vec3f.h>
#include <pxr/base/gf/vec4f.h>
#include <pxr/base/vt/value.h>
#include <pxr/usd/sdf/assetPath.h>
#include <pxr/usd/sdr/shaderNode.h>
#include <pxr/usd/sdr/shaderProperty.h>
#include <pxr/usd/usdShade/material.h>
#include <pxr/usd/usdShade/output.h>
#include <pxr/usd/usdShade/udimUtils.h>
#include <pxr/usd/usdShade/utils.h>
#include <type_traits>

namespace ybi
{

static const pxr::TfToken mtlxToken("mtlx", pxr::TfToken::Immortal);
static const pxr::TfToken mtlxOpenPbr("ND_open_pbr_surface_surfaceshader", pxr::TfToken::Immortal);
static const pxr::TfToken mtlxSurface("ND_surface", pxr::TfToken::Immortal);

namespace mtlxTokens
{
#define MTLX_TOKENS(mtlxName, rendererOpcode) \
    static const pxr::TfToken mtlxName("ND_" #mtlxName, pxr::TfToken::Immortal);

MATERIALX_DEFINITIONS(MTLX_TOKENS)
#undef MTLX_TOKENS

#define MTLX_TOKENS(mtlxName) \
    static const pxr::TfToken mtlxName("ND_" #mtlxName "_bsdf", pxr::TfToken::Immortal);

MATERIALX_BSDFS(MTLX_TOKENS)
#undef MTLX_TOKENS

} // namespace mtlxTokens

static MaterialOpcode ConvertTokenToMaterialOpcode(const pxr::TfToken &token)
{
#define SWITCH(mtlxName, rendererOpcode) \
    if (token == mtlxTokens::mtlxName) \
        return MaterialOpcode::rendererOpcode;

    MATERIALX_DEFINITIONS(SWITCH)
#undef SWITCH

#define SWITCH(mtlxName) \
    if (token == mtlxTokens::mtlxName) \
        return MaterialOpcode::mtlxName##_bsdf;

    MATERIALX_BSDFS(SWITCH)
#undef SWITCH

    YBI_LOGFATAL("MaterialX node %s not supported\n", token.GetText());
    return MaterialOpcode::Invalid;
}

static MaterialOpcode ConvertOutputToMaterialOpcode(const pxr::TfToken &shaderId,
                                                    const pxr::TfToken &outputName)
{
    if (shaderId == pxr::TfToken("ND_artistic_ior"))
    {
        if (outputName == pxr::TfToken("ior"))
            return MaterialOpcode::ArtisticIorIor;
        if (outputName == pxr::TfToken("extinction"))
            return MaterialOpcode::ArtisticIorExtinction;
        YBI_LOGFATAL("Unknown artistic_ior output '%s'\n", outputName.GetText());
        return MaterialOpcode::Invalid;
    }

    return ConvertTokenToMaterialOpcode(shaderId);
}

struct OpenPbrParameterDescriptor
{
    pxr::TfToken token;
    MaterialNodeId OpenPbrMaterial::*member;

    bool (*assignConstant)(MaterialGraph &, MaterialNodeId &, const pxr::VtValue &);
};

template <typename T>
static bool
AssignConstant(MaterialGraph &graph, MaterialNodeId &output, const pxr::VtValue &value);

template <>
bool AssignConstant<float>(MaterialGraph &graph, MaterialNodeId &output, const pxr::VtValue &value)
{
    if (!value.IsHolding<float>())
        return false;

    output = static_cast<MaterialNodeId>(graph.nodes.size());
    graph.nodes.push_back(
        {MaterialOpcode::ConstantFloat, {{"value", value.UncheckedGet<float>()}}});

    return true;
}

template <>
bool AssignConstant<bool>(MaterialGraph &graph, MaterialNodeId &output, const pxr::VtValue &value)
{
    if (!value.IsHolding<bool>())
        return false;

    output = static_cast<MaterialNodeId>(graph.nodes.size());
    graph.nodes.push_back({MaterialOpcode::ConstantBool, {{"value", value.UncheckedGet<bool>()}}});

    return true;
}

template <>
bool AssignConstant<Vec3f>(MaterialGraph &graph, MaterialNodeId &output, const pxr::VtValue &value)
{
    if (!value.IsHolding<pxr::GfVec3f>())
        return false;

    const pxr::GfVec3f &v = value.UncheckedGet<pxr::GfVec3f>();
    output = static_cast<MaterialNodeId>(graph.nodes.size());
    graph.nodes.push_back({MaterialOpcode::ConstantColor, {{"value", Vec3f(v[0], v[1], v[2])}}});

    return true;
}

static const OpenPbrParameterDescriptor openPbrParameterDescriptors[] = {
#define MAKE_PARAMETER_DESCRIPTOR(usdName, type) \
    { \
        pxr::TfToken(#usdName, pxr::TfToken::Immortal), \
        &OpenPbrMaterial::usdName, \
        [](MaterialGraph &graph, MaterialNodeId &output, const pxr::VtValue &value) { \
            return AssignConstant<type>(graph, output, value); \
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
    {
        descriptor.assignConstant(
            material->graph, material->*descriptor.member, property->GetDefaultValueAsSdfType());
    }
}

static void LoadShader(const pxr::UsdShadeShader &shader, pxr::UsdTimeCode timeCode)
{
    const pxr::SdrShaderNodeConstPtr definition = shader.GetShaderNodeForSourceType(mtlxToken);

    if (!definition)
        return;

    for (const pxr::TfToken &inputName : definition->GetShaderInputNames())
    {
        const pxr::UsdShadeInput input = shader.GetInput(inputName);

        if (input)
        {
            const auto producers = input.GetValueProducingAttributes();

            for (const pxr::UsdAttribute &producer : producers)
            {
                printf("help: %s\n", producer.GetName().GetText());
            }

            if (!producers.empty())
            {
                // LoadProducers(producers, timeCode);
                continue;
            }

            pxr::VtValue value;
            if (input.Get(&value, timeCode))
            {
                if (value.IsHolding<int>())
                {
                    int v = value.UncheckedGet<int>();
                    printf("int: %i\n", v);
                }
                continue;
            }
        }

        // The USD prim did not provide this input.
        const pxr::SdrShaderPropertyConstPtr property = definition->GetShaderInput(inputName);

        if (property)
        {
            const pxr::VtValue &defaultValue = property->GetDefaultValueAsSdfType();

            if (!defaultValue.IsEmpty())
            {
            }
        }
    }
}

static void LoadExpression(const pxr::UsdShadeInput &input, pxr::UsdTimeCode timeCode)
{
    auto sources = input.GetConnectedSources();
    for (auto &source : sources)
    {
        pxr::UsdPrim sourcePrim = source.source.GetPrim();
        pxr::UsdShadeShader sourceShader(sourcePrim);

        pxr::TfToken shaderId;
        sourceShader.GetShaderId(&shaderId);

        printf("Node: %s, output: %s\n", shaderId.GetText(), source.sourceName.GetText());
        LoadShader(sourceShader, timeCode);

        for (const pxr::UsdShadeInput &sourceInput : sourceShader.GetInputs())
        {
            LoadExpression(sourceInput, timeCode);
        }
    }
}

static MaterialInputValue LoadConstant(const pxr::UsdAttribute &sourceAttribute,
                                       const pxr::VtValue &value,
                                       pxr::UsdTimeCode timeCode)
{
    if (value.IsHolding<bool>())
    {
        return value.UncheckedGet<bool>();
    }
    else if (value.IsHolding<int>())
    {
        return value.UncheckedGet<int>();
    }
    else if (value.IsHolding<float>())
    {
        return value.UncheckedGet<float>();
    }
    else if (value.IsHolding<pxr::GfVec2f>())
    {
        const pxr::GfVec2f &v = value.UncheckedGet<pxr::GfVec2f>();
        return Vec2f(v[0], v[1]);
    }
    else if (value.IsHolding<pxr::GfVec3f>())
    {
        const pxr::GfVec3f &v = value.UncheckedGet<pxr::GfVec3f>();
        return Vec3f(v[0], v[1], v[2]);
    }
    else if (value.IsHolding<pxr::GfVec4f>())
    {
        const pxr::GfVec4f &v = value.UncheckedGet<pxr::GfVec4f>();
        return Vec4f(v[0], v[1], v[2], v[3]);
    }
    // else if (value.IsHolding<pxr::GfMatrix3d>())
    // {
    //     const double *m = value.UncheckedGet<pxr::GfMatrix3d>().GetArray();
    //     return Matrix3x3{
    //         {{static_cast<float>(m[0]), static_cast<float>(m[1]), static_cast<float>(m[2])},
    //          {static_cast<float>(m[3]), static_cast<float>(m[4]), static_cast<float>(m[5])},
    //          {static_cast<float>(m[6]), static_cast<float>(m[7]), static_cast<float>(m[8])}}};
    // }
    else if (value.IsHolding<pxr::GfMatrix4d>())
    {
        const double *m = value.UncheckedGet<pxr::GfMatrix4d>().GetArray();
        return Matrix4x4(static_cast<float>(m[0]),
                         static_cast<float>(m[1]),
                         static_cast<float>(m[2]),
                         static_cast<float>(m[3]),
                         static_cast<float>(m[4]),
                         static_cast<float>(m[5]),
                         static_cast<float>(m[6]),
                         static_cast<float>(m[7]),
                         static_cast<float>(m[8]),
                         static_cast<float>(m[9]),
                         static_cast<float>(m[10]),
                         static_cast<float>(m[11]),
                         static_cast<float>(m[12]),
                         static_cast<float>(m[13]),
                         static_cast<float>(m[14]),
                         static_cast<float>(m[15]));
    }
    else if (value.IsHolding<std::string>())
    {
        return value.UncheckedGet<std::string>();
    }
    else if (value.IsHolding<pxr::SdfAssetPath>())
    {
        UdimTextureVector filePaths;

        pxr::SdfAssetPath sdfAssetPath = value.UncheckedGet<pxr::SdfAssetPath>();
        ParseUdimPaths(sdfAssetPath, sourceAttribute, timeCode, filePaths);

        return filePaths;
    }
    else if (value.IsHolding<pxr::TfToken>())
    {
        return value.UncheckedGet<pxr::TfToken>().GetString();
    }

    printf("Unsupported MaterialX constant type: %s\n", value.GetTypeName().c_str());
    fflush(stdout);
    YBI_ASSERT(false && "Unsupported MaterialX constant type");
    return false;
}

static void DebugPrintConstant(const MaterialInput &input)
{
    printf("%s: ", input.name.c_str());

    std::visit(
        [](const auto &value) {
            using T = std::decay_t<decltype(value)>;

            if constexpr (std::is_same_v<T, MaterialNodeId>)
            {
                printf("node: %u\n", value);
            }
            else if constexpr (std::is_same_v<T, bool>)
            {
                printf("bool: %s\n", value ? "true" : "false");
            }
            else if constexpr (std::is_same_v<T, int>)
            {
                printf("int: %d\n", value);
            }
            else if constexpr (std::is_same_v<T, float>)
            {
                printf("float: %.9g\n", static_cast<double>(value));
            }
            else if constexpr (std::is_same_v<T, Vec2f>)
            {
                printf("Vec2f: (%.9g, %.9g)\n",
                       static_cast<double>(value.x),
                       static_cast<double>(value.y));
            }
            else if constexpr (std::is_same_v<T, Vec3f>)
            {
                printf("Vec3f: (%.9g, %.9g, %.9g)\n",
                       static_cast<double>(value.x),
                       static_cast<double>(value.y),
                       static_cast<double>(value.z));
            }
            else if constexpr (std::is_same_v<T, Vec4f>)
            {
                printf("Vec4f: (%.9g, %.9g, %.9g, %.9g)\n",
                       static_cast<double>(value.x),
                       static_cast<double>(value.y),
                       static_cast<double>(value.z),
                       static_cast<double>(value.w));
            }
            // else if constexpr (std::is_same_v<T, Matrix3x3>)
            // {
            //     printf("Matrix3x3: ((%.9g, %.9g, %.9g), "
            //            "(%.9g, %.9g, %.9g), "
            //            "(%.9g, %.9g, %.9g))\n",
            //            static_cast<double>(value.m[0][0]),
            //            static_cast<double>(value.m[0][1]),
            //            static_cast<double>(value.m[0][2]),
            //            static_cast<double>(value.m[1][0]),
            //            static_cast<double>(value.m[1][1]),
            //            static_cast<double>(value.m[1][2]),
            //            static_cast<double>(value.m[2][0]),
            //            static_cast<double>(value.m[2][1]),
            //            static_cast<double>(value.m[2][2]));
            // }
            else if constexpr (std::is_same_v<T, Matrix4x4>)
            {
                printf("Matrix4x4: ((%.9g, %.9g, %.9g, %.9g), "
                       "(%.9g, %.9g, %.9g, %.9g), "
                       "(%.9g, %.9g, %.9g, %.9g), "
                       "(%.9g, %.9g, %.9g, %.9g))\n",
                       static_cast<double>(value.m[0][0]),
                       static_cast<double>(value.m[0][1]),
                       static_cast<double>(value.m[0][2]),
                       static_cast<double>(value.m[0][3]),
                       static_cast<double>(value.m[1][0]),
                       static_cast<double>(value.m[1][1]),
                       static_cast<double>(value.m[1][2]),
                       static_cast<double>(value.m[1][3]),
                       static_cast<double>(value.m[2][0]),
                       static_cast<double>(value.m[2][1]),
                       static_cast<double>(value.m[2][2]),
                       static_cast<double>(value.m[2][3]),
                       static_cast<double>(value.m[3][0]),
                       static_cast<double>(value.m[3][1]),
                       static_cast<double>(value.m[3][2]),
                       static_cast<double>(value.m[3][3]));
            }
            else if constexpr (std::is_same_v<T, std::string>)
            {
                printf("string: \"%s\"\n", value.c_str());
            }
            else if constexpr (std::is_same_v<T, pxr::SdfAssetPath>)
            {
                printf("asset: @%s@\n", value.GetAssetPath().c_str());
            }
            else if constexpr (std::is_same_v<T, UdimTextureVector>)
            {
                printf("files: \n");
                for (auto &[filename, udim] : value)
                {
                    printf("%s\n", filename.c_str());
                }
            }
            else
            {
                YBI_ASSERT(0);
            }
        },
        input.value);
}

using NodeMap = std::unordered_map<std::string, MaterialNodeId>;

static MaterialColorSpace GetMaterialColorSpace(const pxr::UsdAttribute &attribute)
{
    const pxr::TfToken colorSpace = attribute.GetColorSpace();
    if (colorSpace.IsEmpty() || colorSpace == pxr::TfToken("raw") ||
        colorSpace == pxr::TfToken("lin_rec709") ||
        colorSpace == pxr::TfToken("lin_rec709_scene") ||
        colorSpace == pxr::TfToken("lin_srgb"))
    {
        return MaterialColorSpace::Linear;
    }
    if (colorSpace == pxr::TfToken("srgb_texture") ||
        colorSpace == pxr::TfToken("srgb_tx") ||
        colorSpace == pxr::TfToken("srgb_rec709_scene"))
    {
        return MaterialColorSpace::SRGB;
    }

    return MaterialColorSpace::Linear;
}

static MaterialNodeId LoadMtlxShaderHelper(const pxr::UsdShadeOutput &shaderOutput,
                                           MaterialGraph &graph,
                                           NodeMap &nodeMap,
                                           pxr::UsdTimeCode timeCode)
{
    const std::string outputPath = shaderOutput.GetAttr().GetPath().GetString();
    if (const auto found = nodeMap.find(outputPath); found != nodeMap.end())
        return found->second;

    pxr::UsdShadeShader shader(shaderOutput.GetPrim());
    pxr::TfToken shaderId;

    if (shader && shader.GetShaderId(&shaderId))
    {
        std::vector<MaterialInput> materialInputs;

        for (const pxr::UsdShadeInput &input : shader.GetInputs())
        {
            const pxr::TfToken inputName = input.GetBaseName();
            printf("input name: %s\n", inputName.GetText());
            const pxr::UsdShadeAttributeVector producers = input.GetValueProducingAttributes();
            YBI_ASSERT(producers.size() <= 1 &&
                       "MaterialX inputs with multiple value producers are unsupported");

            if (producers.empty())
            {
                pxr::VtValue value;
                if (input.Get(&value, timeCode))
                {
                    materialInputs.push_back(
                        {inputName, LoadConstant(input.GetAttr(), value, timeCode)});
                    if (std::holds_alternative<UdimTextureVector>(materialInputs.back().value))
                    {
                        materialInputs.back().colorSpace = GetMaterialColorSpace(input.GetAttr());
                    }
                    DebugPrintConstant(materialInputs.back());
                }
                continue;
            }

            const pxr::UsdAttribute &producer = producers.front();
            const bool isShader = pxr::UsdShadeOutput::IsOutput(producer) &&
                                  producer.GetPrim().IsA<pxr::UsdShadeShader>();

            if (isShader)
            {
                const MaterialNodeId sourceId = LoadMtlxShaderHelper(
                    pxr::UsdShadeOutput(producer), graph, nodeMap, timeCode);
                YBI_ASSERT(sourceId != kInvalidMaterialNodeId);
                materialInputs.push_back({inputName, sourceId});
            }
            else
            {
                pxr::VtValue value;
                if (!producer.Get(&value, timeCode))
                    YBI_LOGFATAL("Unexpected error when handling producer input\n");

                materialInputs.push_back({inputName, LoadConstant(producer, value, timeCode)});
                if (std::holds_alternative<UdimTextureVector>(materialInputs.back().value))
                {
                    const pxr::UsdAttribute &colorSpaceAttribute =
                        input.GetAttr().GetColorSpace().IsEmpty() ? producer : input.GetAttr();
                    materialInputs.back().colorSpace = GetMaterialColorSpace(colorSpaceAttribute);
                }
                DebugPrintConstant(materialInputs.back());
            }
        }

        printf("\n");

        const MaterialNodeId nodeId = static_cast<MaterialNodeId>(graph.nodes.size());

        MaterialOpcode opcode =
            ConvertOutputToMaterialOpcode(shaderId, shaderOutput.GetBaseName());
        YBI_ASSERT(opcode != MaterialOpcode::Invalid);

        graph.nodes.push_back({opcode, std::move(materialInputs)});

        nodeMap[outputPath] = nodeId;
        return nodeId;
    }
    return kInvalidMaterialNodeId;
}

static MaterialNodeId LoadMtlxShader(MaterialGraph &graph,
                                     const pxr::UsdShadeOutput &shadeOutput,
                                     pxr::UsdTimeCode timeCode)
{
    NodeMap nodeMap;
    graph.outputNode = LoadMtlxShaderHelper(shadeOutput, graph, nodeMap, timeCode);

    return graph.outputNode;
}

void LoadOpenPbrMaterial(const pxr::UsdShadeShader &shader,
                         pxr::UsdTimeCode timeCode,
                         OpenPbrMaterial *material)
{
    YBI_ASSERT(shader && material);

    pxr::TfToken idToken;
    YBI_ASSERT(shader.GetShaderId(&idToken) && idToken == mtlxOpenPbr);

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
                descriptor.assignConstant(material->graph, material->*descriptor.member, value);

            continue;
        }

        for (const pxr::UsdAttribute &producer : producers)
        {
            pxr::VtValue value;
            if (pxr::UsdShadeInput::IsInput(producer))
            {
                if (producer.Get(&value, timeCode) &&
                    descriptor.assignConstant(
                        material->graph, material->*descriptor.member, value))
                {
                    break;
                }
                else
                {
                    YBI_LOGFATAL("Unexpected error when handling producer input\n");
                }
            }
            else if (pxr::UsdShadeOutput::IsOutput(producer))
            {
                if (producer.GetPrim().IsA<pxr::UsdShadeShader>())
                {
                    pxr::UsdShadeShader sourceShader(producer.GetPrim());
                    pxr::TfToken shaderId;
                    sourceShader.GetShaderId(&shaderId);

                    material->*descriptor.member =
                        LoadMtlxShader(material->graph, pxr::UsdShadeOutput(producer), timeCode);
                }
                else if (producer.Get(&value, timeCode) &&
                         descriptor.assignConstant(
                             material->graph, material->*descriptor.member, value))
                {
                    break;
                }
                else
                {
                    YBI_LOGFATAL("Unexpected error when handling producer input\n");
                }
            }
        }
    }
}

void LoadMaterialXSurface(const pxr::UsdShadeShader &shader,
                          pxr::UsdTimeCode timeCode,
                          MaterialGraph *material)
{
    YBI_ASSERT(shader && material);
    *material = MaterialGraph{};

    pxr::TfToken idToken;
    YBI_ASSERT(shader.GetShaderId(&idToken) && idToken == mtlxSurface);

    const pxr::UsdShadeInput bsdfInput = shader.GetInput(pxr::TfToken("bsdf"));
    if (!bsdfInput)
        return;

    const pxr::UsdShadeAttributeVector producers = bsdfInput.GetValueProducingAttributes();
    YBI_ASSERT(producers.size() <= 1 &&
               "MaterialX surface inputs with multiple value producers are unsupported");

    for (const pxr::UsdAttribute &producer : producers)
    {
        if (pxr::UsdShadeOutput::IsOutput(producer) &&
            producer.GetPrim().IsA<pxr::UsdShadeShader>())
        {
            pxr::UsdShadeShader sourceShader(producer.GetPrim());
            LoadMtlxShader(*material, pxr::UsdShadeOutput(producer), timeCode);
            return;
        }

        YBI_LOGFATAL("MaterialX surface bsdf input must be connected to a shader output\n");
    }
}

} // namespace ybi
