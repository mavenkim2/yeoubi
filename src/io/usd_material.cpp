#include "io/usd_material.h"
#include "core/array.h"
#include "core/assert.h"
#include "core/matrix4x4.h"
#include "core/vec2f.h"
#include "core/vec4f.h"
#include "scene/mtlx.h"
#include "scene/openpbr_parameters.h"
#include "scene/shader_compilation.h"
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
#include <type_traits>

namespace ybi
{

static const pxr::TfToken mtlxToken("mtlx", pxr::TfToken::Immortal);
static const pxr::TfToken mtlxOpenPbr("ND_open_pbr_surface_surfaceshader", pxr::TfToken::Immortal);

namespace mtlxTokens
{
#define MTLX_TOKENS(name, type) \
    static const pxr::TfToken name("ND_" #name, pxr::TfToken::Immortal);

MATERIALX_DEFINITIONS(MTLX_TOKENS)
#undef MTLX_TOKENS

} // namespace mtlxTokens

static MtlxNodeType ConvertTokenToMtlxNodeType(const pxr::TfToken &token)
{
#define SWITCH(name, type) \
    if (token == mtlxTokens::name) \
        return MtlxNodeType::name;

    MATERIALX_DEFINITIONS(SWITCH)
    YBI_LOGFATAL("MaterialX node %s not supported\n", token.GetText());
    return MtlxNodeType::None;
#undef SWITCH
}

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
                // LoadConstant(inputName, value);
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
                // LoadConstant(inputName, defaultValue);
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

static MaterialInputValue LoadConstant(const pxr::VtValue &value)
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
        return value.UncheckedGet<pxr::SdfAssetPath>().GetAssetPath();
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

            if constexpr (std::is_same_v<T, std::shared_ptr<MaterialNode>>)
            {
                printf("shader: %p\n", static_cast<const void *>(value.get()));
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
        },
        input.value);
}

using NodeMap = std::unordered_map<std::string, std::shared_ptr<MaterialNode>>;

static std::shared_ptr<MaterialNode> LoadMtlxShaderHelper(const pxr::UsdShadeShader &shader,
                                                          NodeMap &nodeMap,
                                                          pxr::UsdTimeCode timeCode)
{
    pxr::TfToken shaderId;

    if (shader && shader.GetShaderId(&shaderId))
    {
        std::string path = shader.GetPath().GetString();
        const pxr::SdrShaderNodeConstPtr definition = shader.GetShaderNodeForSourceType(mtlxToken);

        if (!definition)
            return nullptr;

        std::vector<MaterialInput> materialInputs;

        for (const pxr::TfToken &inputName : definition->GetShaderInputNames())
        {
            const pxr::UsdShadeInput input = shader.GetInput(inputName);
            printf("input name: %s\n", inputName.GetText());
            bool useDefault = true;

            if (input)
            {
                const auto producers = input.GetValueProducingAttributes();
                YBI_ASSERT(producers.size() <= 1 &&
                           "MaterialX inputs with multiple value producers are unsupported");
                pxr::VtValue value;

                if (producers.empty())
                {
                    pxr::VtValue value;
                    if (input.Get(&value, timeCode))
                    {
                        materialInputs.push_back({inputName, LoadConstant(value)});
                        useDefault = false;
                        DebugPrintConstant(materialInputs.back());
                    }
                }

                for (const pxr::UsdAttribute &producer : producers)
                {
                    bool isShader = pxr::UsdShadeOutput::IsOutput(producer) &&
                                    producer.GetPrim().IsA<pxr::UsdShadeShader>();

                    if (isShader)
                    {
                        pxr::UsdShadeShader sourceShader(producer.GetPrim());
                        pxr::TfToken shaderId;
                        sourceShader.GetShaderId(&shaderId);
                        printf("%s\n", shaderId.GetText());

                        if (auto found = nodeMap.find(sourceShader.GetPath().GetString());
                            found == nodeMap.end())
                        {
                            LoadMtlxShaderHelper(sourceShader, nodeMap, timeCode);
                        }

                        auto found = nodeMap.find(sourceShader.GetPath().GetString());
                        YBI_ASSERT(found != nodeMap.end());

                        materialInputs.push_back({inputName, found->second});
                        useDefault = false;
                    }
                    else if (producer.Get(&value, timeCode))
                    {
                        materialInputs.push_back({inputName, LoadConstant(value)});
                        useDefault = false;
                        DebugPrintConstant(materialInputs.back());
                        break;
                    }
                    else
                    {
                        YBI_LOGFATAL("Unexpected error when handling producer input\n");
                    }
                }
            }

            if (useDefault)
            {
                // The USD prim did not provide this input.
                const pxr::SdrShaderPropertyConstPtr property =
                    definition->GetShaderInput(inputName);

                if (property)
                {
                    const pxr::VtValue &defaultValue = property->GetDefaultValueAsSdfType();

                    if (!defaultValue.IsEmpty())
                    {
                        materialInputs.push_back({inputName, LoadConstant(defaultValue)});
                    }
                }
            }
        }

        std::shared_ptr materialNode = std::make_shared<MaterialNode>();
        materialNode->inputs = std::move(materialInputs);
        materialNode->type = ConvertTokenToMtlxNodeType(shaderId);

        nodeMap[path] = materialNode;
        return materialNode;
    }
    return nullptr;
}

static std::shared_ptr<MaterialNode> LoadMtlxShader(const pxr::UsdShadeShader &shader,
                                                    pxr::UsdTimeCode timeCode)
{
    NodeMap nodeMap;
    std::shared_ptr<MaterialNode> result = LoadMtlxShaderHelper(shader, nodeMap, timeCode);
    CompileShader(result);
    printf("\n");

    return result;
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
            if (pxr::UsdShadeInput::IsInput(producer))
            {
                if (producer.Get(&value, timeCode) && descriptor.assignConstant(*material, value))
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

                    printf("input type: %s\n", descriptor.token.GetText());
                    printf("%s\n", shaderId.GetText());
                    LoadMtlxShader(sourceShader, timeCode);
                }
                else if (producer.Get(&value, timeCode) &&
                         descriptor.assignConstant(*material, value))
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

// Shader graph

struct Shader
{
    int inputOffset;
    int outputOffset;
};

} // namespace ybi
