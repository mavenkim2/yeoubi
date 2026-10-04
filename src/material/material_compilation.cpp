#include "material/material_compilation.h"
#include "core/assert.h"
#include "device/device.h"
#include "kernel/bsdf/materialx/bsdf_types.h"
#include "kernel/bsdf/materialx/dielectric_bsdf.h"
#include "material/opcodes.h"
#include <functional>
#include <initializer_list>
#include <iomanip>
#include <limits>
#include <memory>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace ybi
{
namespace
{

static const std::string GET_CLOSURE_PREFIX = "GetClosure_";
static const std::string SET_CLOSURE_PREFIX = "SetClosure_";
static const std::string VALUE_PREFIX = "Value_";
static const std::string INITIALIZE_PREFIX = "Initialize_";
static const std::string SAMPLE_PREFIX = "Sample_";
static const std::string EVALUATE_PREFIX = "Evaluate_";
static const std::string PDF_PREFIX = "Pdf_";
static const std::string ALBEDO_PREFIX = "Albedo_";
static const std::string TRANSMITTANCE_PREFIX = "Transmittance_";

struct MaterialCompilationContext;
struct MaterialNodeBase;

using NodeVector = std::vector<std::unique_ptr<MaterialNodeBase>>;
using SourceNodeMap = std::vector<MaterialNodeBase *>;

using MacroDefinition = std::pair<const char *, std::string>;

static std::string CompileWithMacros(const char *includePath,
                                     std::initializer_list<MacroDefinition> definitions)
{
    std::string result;
    for (const MacroDefinition &definition : definitions)
        result += "#define " + std::string(definition.first) + " " + definition.second + "\n";

    result += "#include \"" + std::string(includePath) + "\"\n";

    for (const MacroDefinition &definition : definitions)
    {
        std::string name = definition.first;
        const size_t parameters = name.find('(');
        if (parameters != std::string::npos)
            name.resize(parameters);
        result += "#undef " + name + "\n";
    }

    return result;
}

static std::string MaterialDebugPrint(const std::string &node,
                                      const char *branch,
                                      MaterialOpcode sampledType,
                                      const char *detailsFormat,
                                      const char *detailsValues)
{
#ifdef YBI_MATERIAL_DEBUG
    return "printf(\"" + node + " " + branch + " -> " +
           ConvertMaterialOpcodeToString(sampledType) +
           " selectedProbability=%f unselectedProbability=%f " + detailsFormat +
           "\\n\", selectedProbability, unselectedProbability, " + detailsValues + ")";
#else
    (void)node;
    (void)branch;
    (void)sampledType;
    (void)detailsFormat;
    (void)detailsValues;
    return "";
#endif
}

static std::string MaterialDebugPrintValues(const std::string &node,
                                            const char *label,
                                            const char *format,
                                            const char *values)
{
#ifdef YBI_MATERIAL_DEBUG
    return "printf(\"" + node + " " + label + " " + format + "\\n\", " + values + ")";
#else
    (void)node;
    (void)label;
    (void)format;
    (void)values;
    return "";
#endif
}

static std::string FloatLiteral(float value)
{
    std::ostringstream stream;
    stream << std::setprecision(std::numeric_limits<float>::max_digits10) << value;

    std::string result = stream.str();
    if (result.find_first_of(".eE") == std::string::npos)
        result += ".0";
    result += "f";
    return result;
}

struct MaterialNodeBase
{
    virtual ~MaterialNodeBase() = default;
    virtual std::string Compile(const MaterialCompilationContext &context,
                                const std::string &materialPrefix) const = 0;
    virtual MaterialOpcode GetNodeType() const = 0;

    virtual bool IsBsdf() const
    {
        return false;
    }

    std::string Name() const
    {
        return "n" + std::to_string(index);
    }

    std::string FullName(const std::string &materialPrefix) const
    {
        return materialPrefix + "_" + Name();
    }

    std::string CompileClosureAccessors(const std::string &materialPrefix,
                                        size_t closureOffset) const
    {
        if (!IsBsdf() && !IsConstantOpcode(GetNodeType()))
        {
            std::string nodeName = FullName(materialPrefix);
            const char *returnType = GetMaterialOpcodeOutputString(GetNodeType());
            return CompileWithMacros(
                "kernel/bsdf/materialx/stubs/closure_impl.h",
                {
                    {"YBI_MATERIAL_VALUE_TYPE", returnType},
                    {"YBI_MATERIAL_SELF_GET_CLOSURE", GET_CLOSURE_PREFIX + nodeName},
                    {"YBI_MATERIAL_SELF_SET_CLOSURE", SET_CLOSURE_PREFIX + nodeName},
                    {"YBI_MATERIAL_VALUE", VALUE_PREFIX + nodeName},
                    {"YBI_MATERIAL_CLOSURE_OFFSET", std::to_string(closureOffset)},
                });
        }
        return "";
    }

    uint32_t index = kInvalidMaterialNodeId;
    std::vector<MaterialNodeBase *> inputs;
};

static std::string ValueCall(const MaterialNodeBase *node, const std::string &materialPrefix)
{
    const bool constant = IsConstantOpcode(node->GetNodeType());
    return (constant ? GET_CLOSURE_PREFIX : VALUE_PREFIX) + node->FullName(materialPrefix) +
           (constant ? "(info)" : "(params, info)");
}

struct BsdfNodeBase : MaterialNodeBase
{
    virtual ~BsdfNodeBase() = default;

    virtual bool IsBsdf() const override
    {
        return true;
    }

    size_t GetClosureSize() const
    {
        size_t count = 0;
        for (MaterialNodeBase *input : inputs)
        {
            if (!input->IsBsdf() && !IsConstantOpcode(input->GetNodeType()))
            {
                count += GetMaterialOpcodeOutputRegisterCount(input->GetNodeType());
            }
        }
        return count;
    }

    std::string CompileInitialization(const MaterialCompilationContext &,
                                      const std::string &materialPrefix) const
    {
        std::string source = "YBI_GPU void " + INITIALIZE_PREFIX + FullName(materialPrefix) +
                             "(const KernelParams &params, "
                             "ClosureInitializationInfo &info)\n{\n    ";

        for (MaterialNodeBase *node : inputs)
        {
            if (node->IsBsdf())
            {
                source += INITIALIZE_PREFIX + node->FullName(materialPrefix) + "(params, info);\n";
            }
            else if (!node->IsBsdf() && !IsConstantOpcode(node->GetNodeType()))
            {
                std::string nodeName = node->FullName(materialPrefix);
                const char *returnType = GetMaterialOpcodeOutputString(node->GetNodeType());

                source += SET_CLOSURE_PREFIX + nodeName + "(params, info);\n";
            }
        }
        source += "}\n";
        return source;
    }

    MaterialNodeBase *normal = nullptr;
};

struct MaterialCompilationContext
{
    explicit MaterialCompilationContext(
        const std::unordered_map<const MaterialNode *, MaterialTextureHandle> &textureNodeToHandle)
        : textureNodeToHandle(textureNodeToHandle)
    {
    }

    const std::unordered_map<const MaterialNode *, MaterialTextureHandle> &textureNodeToHandle;
    NodeVector materialNodes;
    std::vector<MaterialNodeBase *> compilationOrder;

    template <typename T, typename... Args>
    T *AddNode(Args &&...args)
    {
        auto node = std::make_unique<T>(std::forward<Args>(args)...);
        T *result = node.get();
        materialNodes.push_back(std::move(node));
        return result;
    }

    void TopologicalSortHelper(std::unordered_map<MaterialNodeBase *, uint8_t> &state,
                               MaterialNodeBase *node)
    {
        YBI_ASSERT(node);
        uint8_t &nodeState = state[node];
        YBI_ASSERT(nodeState != 1);
        if (nodeState != 0)
            return;

        nodeState = 1;
        for (MaterialNodeBase *input : node->inputs)
            TopologicalSortHelper(state, input);
        nodeState = 2;
        compilationOrder.push_back(node);
    }

    void TopologicalSort(MaterialNodeBase *outputNode)
    {
        compilationOrder.clear();
        compilationOrder.reserve(materialNodes.size());

        std::unordered_map<MaterialNodeBase *, uint8_t> state;
        state.reserve(materialNodes.size());
        TopologicalSortHelper(state, outputNode);

        for (size_t i = 0; i < compilationOrder.size(); ++i)
            compilationOrder[i]->index = static_cast<uint32_t>(compilationOrder.size() - i - 1);
    }
};

struct ConstantFloatNode final : MaterialNodeBase
{
    ConstantFloatNode(float value) : value(value) {}

    std::string Compile(const MaterialCompilationContext &,
                        const std::string &materialPrefix) const override
    {
        return CompileWithMacros(
            "kernel/bsdf/materialx/stubs/constant_impl.h",
            {{"YBI_MATERIAL_VALUE_TYPE", "float"},
             {"YBI_MATERIAL_SELF", GET_CLOSURE_PREFIX + FullName(materialPrefix)},
             {"YBI_MATERIAL_VALUE", FloatLiteral(value)}});
    }

    MaterialOpcode GetNodeType() const override
    {
        return MaterialOpcode::ConstantFloat;
    }

    float value;
};

struct ConstantColor3Node final : MaterialNodeBase
{
    ConstantColor3Node(const Vec3f &value) : value(value) {}

    std::string Compile(const MaterialCompilationContext &,
                        const std::string &materialPrefix) const override
    {
        return CompileWithMacros(
            "kernel/bsdf/materialx/stubs/constant_impl.h",
            {{"YBI_MATERIAL_VALUE_TYPE", "Vec3f"},
             {"YBI_MATERIAL_SELF", GET_CLOSURE_PREFIX + FullName(materialPrefix)},
             {"YBI_MATERIAL_VALUE",
              "Vec3f(" + FloatLiteral(value.x) + ", " + FloatLiteral(value.y) + ", " +
                  FloatLiteral(value.z) + ")"}});
    }

    MaterialOpcode GetNodeType() const override
    {
        return MaterialOpcode::ConstantColor;
    }

    Vec3f value;
};

struct ConstantVector2Node final : MaterialNodeBase
{
    ConstantVector2Node(const Vec2f &value) : value(value) {}

    std::string Compile(const MaterialCompilationContext &,
                        const std::string &materialPrefix) const override
    {
        return CompileWithMacros(
            "kernel/bsdf/materialx/stubs/constant_impl.h",
            {{"YBI_MATERIAL_VALUE_TYPE", "Vec2f"},
             {"YBI_MATERIAL_SELF", GET_CLOSURE_PREFIX + FullName(materialPrefix)},
             {"YBI_MATERIAL_VALUE",
              "Vec2f(" + FloatLiteral(value.x) + ", " + FloatLiteral(value.y) + ")"}});
    }

    MaterialOpcode GetNodeType() const override
    {
        return MaterialOpcode::ConstantVector2;
    }

    Vec2f value;
};

struct ExpressionNode final : MaterialNodeBase
{
    explicit ExpressionNode(MaterialOpcode opcode) : opcode(opcode) {}

    MaterialOpcode GetNodeType() const override
    {
        return opcode;
    }

    std::string Compile(const MaterialCompilationContext &,
                        const std::string &materialPrefix) const override
    {
        const auto arg = [&](size_t index) { return ValueCall(inputs[index], materialPrefix); };
        std::string body;
        switch (opcode)
        {
            case MaterialOpcode::ViewDirectionVector3:
                body = "return -info.rayDirection;";
                break;
            case MaterialOpcode::NormalVector3:
                body = "return info.shadingNormal;";
                break;
            case MaterialOpcode::DotProductVector3:
                body = "const float value = Dot(" + arg(0) + ", " + arg(1) + "); ";
#ifdef YBI_MATERIAL_DEBUG
                if ((inputs[0]->GetNodeType() == MaterialOpcode::ViewDirectionVector3 &&
                     inputs[1]->GetNodeType() == MaterialOpcode::NormalVector3) ||
                    (inputs[0]->GetNodeType() == MaterialOpcode::NormalVector3 &&
                     inputs[1]->GetNodeType() == MaterialOpcode::ViewDirectionVector3))
                {
                    body += MaterialDebugPrintValues(
                                FullName(materialPrefix), "NdotV", "value=%f", "value") +
                            "; ";
                }
#endif
                body += "return value;";
                break;
            case MaterialOpcode::AbsValFloat:
                body = "return Abs(" + arg(0) + ");";
                break;
            case MaterialOpcode::AsinFloat:
                body = "return asinf(" + arg(0) + ");";
                break;
            case MaterialOpcode::PowerFloat:
                body = "return powf(" + arg(0) + ", " + arg(1) + ");";
                break;
            case MaterialOpcode::DivideFloat:
            case MaterialOpcode::DivideVector3:
                body = "return " + arg(0) + " / " + arg(1) + ";";
                break;
            case MaterialOpcode::MultiplyFloat:
            case MaterialOpcode::MultiplyVector3:
                body = "return " + arg(0) + " * " + arg(1) + ";";
                break;
            case MaterialOpcode::ExpVector3:
                body = "const Vec3f value = " + arg(0) +
                       "; return Vec3f(expf(value.x), expf(value.y), expf(value.z));";
                break;
            case MaterialOpcode::LnVector3:
                body = "const Vec3f value = " + arg(0) +
                       "; return Vec3f(logf(value.x), logf(value.y), logf(value.z));";
                break;
            case MaterialOpcode::NormalizeVector3:
                body = "const Vec3f value = " + arg(0) +
                       "; return LengthSquared(value) > 0.f ? Normalize(value) : Vec3f(0.f);";
                break;
            case MaterialOpcode::ConvertFloatVector2:
                body = "return Vec2f(" + arg(0) + ");";
                break;
            case MaterialOpcode::ConvertFloatVector3:
                body = "return Vec3f(" + arg(0) + ");";
                break;
            case MaterialOpcode::ConvertVector3Color3:
            case MaterialOpcode::ConvertColor3Vector3:
                body = "return " + arg(0) + ";";
                break;
            case MaterialOpcode::ExtractVector3:
                body = "return " + arg(0) + "[" + std::to_string(component) + "];";
                break;
            case MaterialOpcode::AddFloat:
                body = "return " + arg(0) + " + " + arg(1) + ";";
                break;
            case MaterialOpcode::AddVector3:
#ifdef YBI_MATERIAL_DEBUG
                if (inputs[0]->GetNodeType() == MaterialOpcode::NormalMapFloat &&
                    inputs[1]->GetNodeType() == MaterialOpcode::NormalMapFloat)
                {
                    body =
                        "const Vec3f baseNormal = " + arg(0) +
                        "; const Vec3f addedNormal = " + arg(1) + "; " +
                        MaterialDebugPrintValues(FullName(materialPrefix),
                                                 "add_vector3 normalmap in2 beforeAdd",
                                                 "normal=(%f,%f,%f)",
                                                 "addedNormal.x, addedNormal.y, addedNormal.z") +
                        "; return baseNormal + addedNormal;";
                    break;
                }
#endif
                body = "return " + arg(0) + " + " + arg(1) + ";";
                break;
            case MaterialOpcode::InvertFloat:
                body = "return " + arg(1) + " - " + arg(0) + ";";
                break;
            case MaterialOpcode::MixFloat:
            case MaterialOpcode::MixVector2:
            case MaterialOpcode::MixVector3:
                body = "return Lerp(" + arg(0) + ", " + arg(1) + ", " + arg(2) + ");";
                break;
            case MaterialOpcode::ContrastFloat:
                body = "const float value = " + arg(0) + "; const float amount = " + arg(1) +
                       "; const float pivot = " + arg(2) +
                       "; return (value - pivot) * amount + pivot;";
                break;
            case MaterialOpcode::ClampFloat:
                body = "return Clamp(" + arg(0) + ", " + arg(1) + ", " + arg(2) + ");";
                break;
            default:
                YBI_ASSERT(false && "Unsupported MaterialX value expression");
                return {};
        }

        return CompileWithMacros(
            "kernel/bsdf/materialx/stubs/value_impl.h",
            {{"YBI_MATERIAL_VALUE_TYPE", GetMaterialOpcodeOutputString(opcode)},
             {"YBI_MATERIAL_SELF", VALUE_PREFIX + FullName(materialPrefix)},
             {"YBI_MATERIAL_VALUE_BODY", body}});
    }

    MaterialOpcode opcode;
    int component = 0;
};

struct RemapNode final : MaterialNodeBase
{
    RemapNode(MaterialValueType valueType) : valueType(valueType) {}

    std::string Compile(const MaterialCompilationContext &context,
                        const std::string &materialPrefix) const override
    {
        return CompileWithMacros(
            "kernel/bsdf/materialx/stubs/remap_impl.h",
            {
                {"YBI_VALUE", ValueCall(value, materialPrefix)},
                {"YBI_IN_LOW", ValueCall(inLow, materialPrefix)},
                {"YBI_IN_HIGH", ValueCall(inHigh, materialPrefix)},
                {"YBI_OUT_LOW", ValueCall(outLow, materialPrefix)},
                {"YBI_OUT_HIGH", ValueCall(outHigh, materialPrefix)},
                {"YBI_SELF", VALUE_PREFIX + FullName(materialPrefix)},
                {"YBI_VALUE_TYPE", std::string(GetMaterialValueTypeString(valueType))},
            });
    }

    MaterialOpcode GetNodeType() const override
    {
        switch (valueType)
        {
            case MaterialValueType::Vec3:
            {
                return MaterialOpcode::RemapVector3;
            }
            case MaterialValueType::Float:
                return MaterialOpcode::RemapFloat;
            default:
            {
                YBI_ASSERT(0);
                return MaterialOpcode::Invalid;
            }
        }
    }

    MaterialNodeBase *value;
    MaterialNodeBase *inLow;
    MaterialNodeBase *inHigh;
    MaterialNodeBase *outLow;
    MaterialNodeBase *outHigh;

    MaterialValueType valueType;
};

struct NormalMapNode final : MaterialNodeBase
{
    MaterialOpcode GetNodeType() const override
    {
        return MaterialOpcode::NormalMapFloat;
    }

    std::string Compile(const MaterialCompilationContext &,
                        const std::string &materialPrefix) const override
    {
        return CompileWithMacros(
            "kernel/bsdf/materialx/stubs/normalmap_impl.h",
            {{"YBI_MATERIAL_SELF", VALUE_PREFIX + FullName(materialPrefix)},
             {"YBI_MATERIAL_INPUT", ValueCall(value, materialPrefix)},
             {"YBI_MATERIAL_SCALE", ValueCall(scale, materialPrefix)},
             {"YBI_MATERIAL_NORMAL", ValueCall(normal, materialPrefix)},
             {"YBI_MATERIAL_TANGENT", tangent ? ValueCall(tangent, materialPrefix) : "info.dpdu"},
             {"YBI_MATERIAL_BITANGENT",
              bitangent ? ValueCall(bitangent, materialPrefix) : "Cross(normal, tangent)"}});
    }

    MaterialNodeBase *value = nullptr;
    MaterialNodeBase *scale = nullptr;
    MaterialNodeBase *normal = nullptr;
    MaterialNodeBase *tangent = nullptr;
    MaterialNodeBase *bitangent = nullptr;
};

struct ArtisticIorIorNode final : MaterialNodeBase
{
    MaterialOpcode GetNodeType() const override
    {
        return MaterialOpcode::ArtisticIorIor;
    }

    std::string Compile(const MaterialCompilationContext &,
                        const std::string &materialPrefix) const override
    {
        return CompileWithMacros(
            "kernel/bsdf/materialx/stubs/artistic_ior_ior_impl.h",
            {{"YBI_MATERIAL_SELF", VALUE_PREFIX + FullName(materialPrefix)},
             {"YBI_MATERIAL_REFLECTIVITY", ValueCall(reflectivity, materialPrefix)},
             {"YBI_MATERIAL_EDGE_COLOR", ValueCall(edgeColor, materialPrefix)},
             {"YBI_MATERIAL_DEBUG_PRINT(inputReflectivity, edgeColor, output)",
              MaterialDebugPrintValues(
                  FullName(materialPrefix),
                  "artistic_ior ior",
                  "reflectivity=(%f,%f,%f) edgeColor=(%f,%f,%f) output=(%f,%f,%f)",
                  "inputReflectivity.x, inputReflectivity.y, inputReflectivity.z, "
                  "edgeColor.x, edgeColor.y, edgeColor.z, output.x, output.y, output.z")}});
    }

    MaterialNodeBase *reflectivity = nullptr;
    MaterialNodeBase *edgeColor = nullptr;
};

struct ArtisticIorExtinctionNode final : MaterialNodeBase
{
    MaterialOpcode GetNodeType() const override
    {
        return MaterialOpcode::ArtisticIorExtinction;
    }

    std::string Compile(const MaterialCompilationContext &,
                        const std::string &materialPrefix) const override
    {
        return CompileWithMacros(
            "kernel/bsdf/materialx/stubs/artistic_ior_extinction_impl.h",
            {{"YBI_MATERIAL_SELF", VALUE_PREFIX + FullName(materialPrefix)},
             {"YBI_MATERIAL_REFLECTIVITY", ValueCall(reflectivity, materialPrefix)},
             {"YBI_MATERIAL_EDGE_COLOR", ValueCall(edgeColor, materialPrefix)},
             {"YBI_MATERIAL_DEBUG_PRINT(inputReflectivity, edgeColor, output)",
              MaterialDebugPrintValues(
                  FullName(materialPrefix),
                  "artistic_ior extinction",
                  "reflectivity=(%f,%f,%f) edgeColor=(%f,%f,%f) output=(%f,%f,%f)",
                  "inputReflectivity.x, inputReflectivity.y, inputReflectivity.z, "
                  "edgeColor.x, edgeColor.y, edgeColor.z, output.x, output.y, output.z")}});
    }

    MaterialNodeBase *reflectivity = nullptr;
    MaterialNodeBase *edgeColor = nullptr;
};

struct ImageNode final : MaterialNodeBase
{
    ImageNode(MaterialValueType valueType, const MaterialNode *sourceNode)
        : valueType(valueType), sourceNode(sourceNode)
    {
    }

    std::string Compile(const MaterialCompilationContext &context,
                        const std::string &materialPrefix) const override
    {
        const auto found = context.textureNodeToHandle.find(sourceNode);
        YBI_ASSERT(found != context.textureNodeToHandle.end());
        const MaterialTextureHandle handle = found->second;

        const char *returnType = nullptr;
        const char *returnExpression = nullptr;
        switch (valueType)
        {
            case MaterialValueType::Float:
                returnType = "float";
                returnExpression = "value.x";
                break;
            case MaterialValueType::Vec2:
                returnType = "Vec2f";
                returnExpression = "value.Xy()";
                break;
            case MaterialValueType::Vec3:
                returnType = "Vec3f";
                returnExpression = "value.Xyz()";
                break;
            case MaterialValueType::Vec4:
                returnType = "Vec4f";
                returnExpression = "value";
                break;
            default:
                YBI_ASSERT(false && "Unsupported MaterialX image output type");
                return {};
        }

        const std::string handleLiteral = "MaterialTextureHandle{" +
                                          std::to_string(handle.udimBase) + "u, " +
                                          std::to_string(handle.udimTableOffset) + "u}";

        const bool tiled = sourceNode->type == MaterialOpcode::TiledImageFloat ||
                           sourceNode->type == MaterialOpcode::TiledImageVector2 ||
                           sourceNode->type == MaterialOpcode::TiledImageVector3 ||
                           sourceNode->type == MaterialOpcode::TiledImageVector4;
        if (tiled)
        {
            return CompileWithMacros(
                "kernel/bsdf/materialx/stubs/tiledimage_impl.h",
                {{"YBI_MATERIAL_RETURN", returnType},
                 {"YBI_MATERIAL_SELF", VALUE_PREFIX + FullName(materialPrefix)},
                 {"YBI_MATERIAL_HANDLE", handleLiteral},
                 {"YBI_MATERIAL_RESULT", returnExpression},
                 {"YBI_MATERIAL_UV_TILING", ValueCall(uvTiling, materialPrefix)},
                 {"YBI_MATERIAL_UV_OFFSET", ValueCall(uvOffset, materialPrefix)},
                 {"YBI_MATERIAL_REALWORLD_IMAGE_SIZE",
                  ValueCall(realWorldImageSize, materialPrefix)},
                 {"YBI_MATERIAL_REALWORLD_TILE_SIZE",
                  ValueCall(realWorldTileSize, materialPrefix)}});
        }

        return CompileWithMacros("kernel/bsdf/materialx/stubs/image_impl.h",
                                 {{"YBI_MATERIAL_RETURN", returnType},
                                  {"YBI_MATERIAL_SELF", VALUE_PREFIX + FullName(materialPrefix)},
                                  {"YBI_MATERIAL_HANDLE", handleLiteral},
                                  {"YBI_MATERIAL_RESULT", returnExpression}});
    }

    MaterialOpcode GetNodeType() const override
    {
        return sourceNode->type;
    }

    MaterialValueType valueType;
    // Temporary bridge to the texture handle assigned before source compilation.
    const MaterialNode *sourceNode;
    MaterialNodeBase *uvTiling = nullptr;
    MaterialNodeBase *uvOffset = nullptr;
    MaterialNodeBase *realWorldImageSize = nullptr;
    MaterialNodeBase *realWorldTileSize = nullptr;
};

struct SheenBsdfNode final : BsdfNodeBase
{
    std::string Compile(const MaterialCompilationContext &,
                        const std::string &materialPrefix) const override
    {
        const std::string self = FullName(materialPrefix);
        const std::string weightName = weight->FullName(materialPrefix);
        const std::string colorName = color->FullName(materialPrefix);
        const std::string roughnessName = roughness->FullName(materialPrefix);

        return CompileWithMacros(
            "kernel/bsdf/materialx/stubs/sheen_impl.h",
            {{"YBI_MATERIAL_SELF_SAMPLE", SAMPLE_PREFIX + self},
             {"YBI_MATERIAL_SELF_EVALUATE", EVALUATE_PREFIX + self},
             {"YBI_MATERIAL_SELF_PDF", PDF_PREFIX + self},
             {"YBI_MATERIAL_SELF_ALBEDO", ALBEDO_PREFIX + self},
             {"YBI_MATERIAL_SELF_TRANSMITTANCE", TRANSMITTANCE_PREFIX + self},
             {"YBI_MATERIAL_DEBUG_PRINT_TRANSMITTANCE(weight, albedo, transmittance)",
              MaterialDebugPrintValues(self,
                                       "sheen transmittance",
                                       "weight=%f albedo=%f value=%f",
                                       "weight, albedo, transmittance")},
             {"YBI_MATERIAL_WEIGHT", GET_CLOSURE_PREFIX + weightName},
             {"YBI_MATERIAL_COLOR", GET_CLOSURE_PREFIX + colorName},
             {"YBI_MATERIAL_ROUGHNESS", GET_CLOSURE_PREFIX + roughnessName},
             {"YBI_MATERIAL_NORMAL", GET_CLOSURE_PREFIX + normal->FullName(materialPrefix)}});
    }

    MaterialOpcode GetNodeType() const override
    {
        return MaterialOpcode::sheen_bsdf;
    }

    MaterialNodeBase *weight = nullptr;
    MaterialNodeBase *color = nullptr;
    MaterialNodeBase *roughness = nullptr;
};

struct DielectricBsdfNode final : BsdfNodeBase
{
    std::string Compile(const MaterialCompilationContext &,
                        const std::string &materialPrefix) const override
    {
        const std::string self = FullName(materialPrefix);
        const std::string weightName = weight->FullName(materialPrefix);
        const std::string tintName = tint->FullName(materialPrefix);
        const std::string iorName = ior->FullName(materialPrefix);
        const std::string roughnessName = roughness->FullName(materialPrefix);

        return CompileWithMacros(
            "kernel/bsdf/materialx/stubs/dielectric_impl.h",
            {{"YBI_MATERIAL_SELF_SAMPLE", SAMPLE_PREFIX + self},
             {"YBI_MATERIAL_SELF_EVALUATE", EVALUATE_PREFIX + self},
             {"YBI_MATERIAL_SELF_PDF", PDF_PREFIX + self},
             {"YBI_MATERIAL_SELF_ALBEDO", ALBEDO_PREFIX + self},
             {"YBI_MATERIAL_SELF_TRANSMITTANCE", TRANSMITTANCE_PREFIX + self},
             {"YBI_MATERIAL_DEBUG_PRINT_TRANSMITTANCE(weight, albedo, transmittance)",
              MaterialDebugPrintValues(self,
                                       "dielectric transmittance",
                                       "weight=%f albedo=%f value=%f",
                                       "weight, albedo, transmittance")},
             {"YBI_MATERIAL_WEIGHT", GET_CLOSURE_PREFIX + weightName},
             {"YBI_MATERIAL_TINT", GET_CLOSURE_PREFIX + tintName},
             {"YBI_MATERIAL_IOR", GET_CLOSURE_PREFIX + iorName},
             {"YBI_MATERIAL_ROUGHNESS", GET_CLOSURE_PREFIX + roughnessName},
             {"YBI_MATERIAL_NORMAL", GET_CLOSURE_PREFIX + normal->FullName(materialPrefix)}});
    }

    MaterialOpcode GetNodeType() const override
    {
        return MaterialOpcode::dielectric_bsdf;
    }

    MaterialNodeBase *weight = nullptr;
    MaterialNodeBase *tint = nullptr;
    MaterialNodeBase *ior = nullptr;
    MaterialNodeBase *roughness = nullptr;
};

struct OrenNayarBsdfNode final : BsdfNodeBase
{
    std::string Compile(const MaterialCompilationContext &,
                        const std::string &materialPrefix) const override
    {
        const std::string self = FullName(materialPrefix);
        const std::string weightName = weight->FullName(materialPrefix);
        const std::string colorName = color->FullName(materialPrefix);
        const std::string roughnessName = roughness->FullName(materialPrefix);

        return CompileWithMacros(
            "kernel/bsdf/materialx/stubs/oren_nayar_impl.h",
            {{"YBI_MATERIAL_SELF_SAMPLE", SAMPLE_PREFIX + self},
             {"YBI_MATERIAL_SELF_EVALUATE", EVALUATE_PREFIX + self},
             {"YBI_MATERIAL_SELF_PDF", PDF_PREFIX + self},
             {"YBI_MATERIAL_SELF_ALBEDO", ALBEDO_PREFIX + self},
             {"YBI_MATERIAL_SELF_TRANSMITTANCE", TRANSMITTANCE_PREFIX + self},
             {"YBI_MATERIAL_DEBUG_PRINT_TRANSMITTANCE(weight, transmittance)",
              MaterialDebugPrintValues(self,
                                       "oren_nayar transmittance",
                                       "weight=%f value=%f",
                                       "weight, transmittance")},
             {"YBI_MATERIAL_WEIGHT", GET_CLOSURE_PREFIX + weightName},
             {"YBI_MATERIAL_COLOR", GET_CLOSURE_PREFIX + colorName},
             {"YBI_MATERIAL_ROUGHNESS", GET_CLOSURE_PREFIX + roughnessName},
             {"YBI_MATERIAL_NORMAL", GET_CLOSURE_PREFIX + normal->FullName(materialPrefix)}});
    }

    MaterialOpcode GetNodeType() const override
    {
        return MaterialOpcode::oren_nayar_diffuse_bsdf;
    }

    MaterialNodeBase *weight = nullptr;
    MaterialNodeBase *color = nullptr;
    MaterialNodeBase *roughness = nullptr;
};

struct ConductorBsdfNode final : BsdfNodeBase
{
    std::string Compile(const MaterialCompilationContext &,
                        const std::string &materialPrefix) const override
    {
        const std::string self = FullName(materialPrefix);
        return CompileWithMacros(
            "kernel/bsdf/materialx/stubs/conductor_impl.h",
            {{"YBI_MATERIAL_SELF_SAMPLE", SAMPLE_PREFIX + self},
             {"YBI_MATERIAL_SELF_EVALUATE", EVALUATE_PREFIX + self},
             {"YBI_MATERIAL_SELF_PDF", PDF_PREFIX + self},
             {"YBI_MATERIAL_SELF_ALBEDO", ALBEDO_PREFIX + self},
             {"YBI_MATERIAL_SELF_TRANSMITTANCE", TRANSMITTANCE_PREFIX + self},
             {"YBI_MATERIAL_WEIGHT", GET_CLOSURE_PREFIX + weight->FullName(materialPrefix)},
             {"YBI_MATERIAL_IOR", GET_CLOSURE_PREFIX + ior->FullName(materialPrefix)},
             {"YBI_MATERIAL_EXTINCTION",
              GET_CLOSURE_PREFIX + extinction->FullName(materialPrefix)},
             {"YBI_MATERIAL_ROUGHNESS", GET_CLOSURE_PREFIX + roughness->FullName(materialPrefix)},
             {"YBI_MATERIAL_NORMAL", GET_CLOSURE_PREFIX + normal->FullName(materialPrefix)}});
    }

    MaterialOpcode GetNodeType() const override
    {
        return MaterialOpcode::conductor_bsdf;
    }

    MaterialNodeBase *weight = nullptr;
    MaterialNodeBase *ior = nullptr;
    MaterialNodeBase *extinction = nullptr;
    MaterialNodeBase *roughness = nullptr;
};

struct MixBsdfNode final : BsdfNodeBase
{
    std::string Compile(const MaterialCompilationContext &,
                        const std::string &materialPrefix) const override
    {
        const std::string self = FullName(materialPrefix);
        const std::string bgName = bg->FullName(materialPrefix);
        const std::string fgName = fg->FullName(materialPrefix);
        return CompileWithMacros(
            "kernel/bsdf/materialx/stubs/mix_impl.h",
            {{"YBI_MATERIAL_SELF_SAMPLE", SAMPLE_PREFIX + self},
             {"YBI_MATERIAL_SELF_EVALUATE", EVALUATE_PREFIX + self},
             {"YBI_MATERIAL_SELF_PDF", PDF_PREFIX + self},
             {"YBI_MATERIAL_SELF_ALBEDO", ALBEDO_PREFIX + self},
             {"YBI_MATERIAL_SELF_TRANSMITTANCE", TRANSMITTANCE_PREFIX + self},
             {"YBI_MATERIAL_DEBUG_PRINT_TRANSMITTANCE(mix, bgTransmittance, fgTransmittance, "
              "transmittance)",
              MaterialDebugPrintValues(self,
                                       "mix transmittance",
                                       "mix=%f bg=%f fg=%f value=%f",
                                       "mix, bgTransmittance, fgTransmittance, transmittance")},
             {"YBI_MATERIAL_MIX", GET_CLOSURE_PREFIX + mix->FullName(materialPrefix)},
             {"YBI_MATERIAL_DEBUG_PRINT_0(selectedProbability, unselectedProbability, bgWeight, "
              "fgWeight)",
              MaterialDebugPrint(self,
                                 "mix bg",
                                 bg->GetNodeType(),
                                 "bgWeight=%f fgWeight=%f",
                                 "bgWeight, fgWeight")},
             {"YBI_MATERIAL_DEBUG_PRINT_1(selectedProbability, unselectedProbability, bgWeight, "
              "fgWeight)",
              MaterialDebugPrint(self,
                                 "mix fg",
                                 fg->GetNodeType(),
                                 "bgWeight=%f fgWeight=%f",
                                 "bgWeight, fgWeight")},
             {"YBI_MATERIAL_SAMPLE_BG", SAMPLE_PREFIX + bgName},
             {"YBI_MATERIAL_EVALUATE_BG", EVALUATE_PREFIX + bgName},
             {"YBI_MATERIAL_PDF_BG", PDF_PREFIX + bgName},
             {"YBI_MATERIAL_ALBEDO_BG", ALBEDO_PREFIX + bgName},
             {"YBI_MATERIAL_TRANSMITTANCE_BG", TRANSMITTANCE_PREFIX + bgName},
             {"YBI_MATERIAL_SAMPLE_FG", SAMPLE_PREFIX + fgName},
             {"YBI_MATERIAL_EVALUATE_FG", EVALUATE_PREFIX + fgName},
             {"YBI_MATERIAL_PDF_FG", PDF_PREFIX + fgName},
             {"YBI_MATERIAL_ALBEDO_FG", ALBEDO_PREFIX + fgName},
             {"YBI_MATERIAL_TRANSMITTANCE_FG", TRANSMITTANCE_PREFIX + fgName}});
    }

    MaterialOpcode GetNodeType() const override
    {
        return MaterialOpcode::mix_bsdf;
    }

    MaterialNodeBase *bg = nullptr;
    MaterialNodeBase *fg = nullptr;
    MaterialNodeBase *mix = nullptr;
};

struct LayerBsdfNode final : BsdfNodeBase
{
    std::string Compile(const MaterialCompilationContext &,
                        const std::string &materialPrefix) const override
    {
        const std::string self = FullName(materialPrefix);
        const std::string topName = top->FullName(materialPrefix);
        const std::string baseName = base->FullName(materialPrefix);

        return CompileWithMacros(
            "kernel/bsdf/materialx/stubs/layer_impl.h",
            {{"YBI_MATERIAL_SELF_SAMPLE", SAMPLE_PREFIX + self},
             {"YBI_MATERIAL_SELF_EVALUATE", EVALUATE_PREFIX + self},
             {"YBI_MATERIAL_SELF_PDF", PDF_PREFIX + self},
             {"YBI_MATERIAL_SELF_ALBEDO", ALBEDO_PREFIX + self},
             {"YBI_MATERIAL_SELF_TRANSMITTANCE", TRANSMITTANCE_PREFIX + self},
             {"YBI_MATERIAL_DEBUG_PRINT_TRANSMITTANCE(topTransmittance, bottomTransmittance, "
              "transmittance)",
              MaterialDebugPrintValues(self,
                                       "layer transmittance",
                                       "top=%f bottom=%f value=%f",
                                       "topTransmittance, bottomTransmittance, transmittance")},
             {"YBI_MATERIAL_DEBUG_PRINT_0(selectedProbability, unselectedProbability, topAlbedo, "
              "transmittanceTop, bottomAlbedo, weightedBottomAlbedo)",
              MaterialDebugPrint(
                  self,
                  "layer top",
                  top->GetNodeType(),
                  "topAlbedo=%f transmittanceTop=%f bottomAlbedo=%f weightedBottomAlbedo=%f",
                  "topAlbedo, transmittanceTop, bottomAlbedo, weightedBottomAlbedo")},
             {"YBI_MATERIAL_DEBUG_PRINT_1(selectedProbability, unselectedProbability, topAlbedo, "
              "transmittanceTop, bottomAlbedo, weightedBottomAlbedo)",
              MaterialDebugPrint(
                  self,
                  "layer bottom",
                  base->GetNodeType(),
                  "topAlbedo=%f transmittanceTop=%f bottomAlbedo=%f weightedBottomAlbedo=%f",
                  "topAlbedo, transmittanceTop, bottomAlbedo, weightedBottomAlbedo")},
             {"YBI_MATERIAL_SAMPLE_TOP", SAMPLE_PREFIX + topName},
             {"YBI_MATERIAL_EVALUATE_TOP", EVALUATE_PREFIX + topName},
             {"YBI_MATERIAL_PDF_TOP", PDF_PREFIX + topName},
             {"YBI_MATERIAL_ALBEDO_TOP", ALBEDO_PREFIX + topName},
             {"YBI_MATERIAL_TRANSMITTANCE_TOP", TRANSMITTANCE_PREFIX + topName},
             {"YBI_MATERIAL_SAMPLE_BOTTOM", SAMPLE_PREFIX + baseName},
             {"YBI_MATERIAL_EVALUATE_BOTTOM", EVALUATE_PREFIX + baseName},
             {"YBI_MATERIAL_PDF_BOTTOM", PDF_PREFIX + baseName},
             {"YBI_MATERIAL_ALBEDO_BOTTOM", ALBEDO_PREFIX + baseName},
             {"YBI_MATERIAL_TRANSMITTANCE_BOTTOM", TRANSMITTANCE_PREFIX + baseName}});
    }

    MaterialOpcode GetNodeType() const override
    {
        return MaterialOpcode::layer_bsdf;
    }

    MaterialNodeBase *top = nullptr;
    MaterialNodeBase *base = nullptr;
};

static const MaterialInput *FindInput(const MaterialNode &node, const char *name)
{
    for (const MaterialInput &input : node.inputs)
    {
        if (input.name == name)
            return &input;
    }
    return nullptr;
}

static MaterialNodeBase *GetMappedNode(const SourceNodeMap &sourceNodes, MaterialNodeId node)
{
    YBI_ASSERT(node < sourceNodes.size());
    YBI_ASSERT(sourceNodes[node]);
    return sourceNodes[node];
}

static MaterialNodeBase *
GetNodeInput(const SourceNodeMap &sourceNodes, const MaterialNode &node, const char *name)
{
    const MaterialInput *input = FindInput(node, name);
    YBI_ASSERT(input);
    const MaterialNodeId *connectedNode = std::get_if<MaterialNodeId>(&input->value);
    YBI_ASSERT(connectedNode);
    return GetMappedNode(sourceNodes, *connectedNode);
}

static MaterialNodeBase *AddFloatInput(MaterialCompilationContext &context,
                                       const SourceNodeMap &sourceNodes,
                                       const MaterialNode &node,
                                       const char *name,
                                       float defaultValue)
{
    const MaterialInput *input = FindInput(node, name);
    if (!input)
        return context.AddNode<ConstantFloatNode>(defaultValue);

    if (const MaterialNodeId *connectedNode = std::get_if<MaterialNodeId>(&input->value))
        return GetMappedNode(sourceNodes, *connectedNode);

    const float *value = std::get_if<float>(&input->value);
    YBI_ASSERT(value);
    return context.AddNode<ConstantFloatNode>(*value);
}

static MaterialNodeBase *AddColor3Input(MaterialCompilationContext &context,
                                        const SourceNodeMap &sourceNodes,
                                        const MaterialNode &node,
                                        const char *name,
                                        const Vec3f &defaultValue)
{
    const MaterialInput *input = FindInput(node, name);
    if (!input)
        return context.AddNode<ConstantColor3Node>(defaultValue);

    if (const MaterialNodeId *connectedNode = std::get_if<MaterialNodeId>(&input->value))
        return GetMappedNode(sourceNodes, *connectedNode);

    const Vec3f *value = std::get_if<Vec3f>(&input->value);
    YBI_ASSERT(value);
    return context.AddNode<ConstantColor3Node>(*value);
}

static MaterialNodeBase *AddNormalInput(MaterialCompilationContext &context,
                                        const SourceNodeMap &sourceNodes,
                                        const MaterialNode &node,
                                        MaterialNodeBase *defaultNormal)
{
    if (!FindInput(node, "normal"))
        return defaultNormal;
    return AddColor3Input(context, sourceNodes, node, "normal", Vec3f(0.f));
}

static MaterialNodeBase *AddVector2Input(MaterialCompilationContext &context,
                                         const SourceNodeMap &sourceNodes,
                                         const MaterialNode &node,
                                         const char *name,
                                         const Vec2f &defaultValue)
{
    const MaterialInput *input = FindInput(node, name);
    if (!input)
        return context.AddNode<ConstantVector2Node>(defaultValue);

    if (const MaterialNodeId *connectedNode = std::get_if<MaterialNodeId>(&input->value))
        return GetMappedNode(sourceNodes, *connectedNode);

    const Vec2f *value = std::get_if<Vec2f>(&input->value);
    YBI_ASSERT(value);
    return context.AddNode<ConstantVector2Node>(*value);
}

static void
AssertStringInput(const MaterialNode &node, const char *name, const char *expectedValue)
{
    const MaterialInput *input = FindInput(node, name);
    if (!input)
        return;

    const std::string *value = std::get_if<std::string>(&input->value);
    YBI_ASSERT(value && *value == expectedValue);
}

static void AssertZeroFloatInput(const MaterialNode &node, const char *name)
{
    const MaterialInput *input = FindInput(node, name);
    if (!input)
        return;

    const float *value = std::get_if<float>(&input->value);
    YBI_ASSERT(value && *value == 0.f);
}

static void AssertFalseBoolInput(const MaterialNode &node, const char *name)
{
    const MaterialInput *input = FindInput(node, name);
    if (!input)
        return;

    const bool *value = std::get_if<bool>(&input->value);
    YBI_ASSERT(value && !*value);
}

static std::vector<bool> FindReachableNodes(const MaterialGraph &graph)
{
    std::vector<bool> reachable(graph.nodes.size(), false);

    std::function<void(MaterialNodeId)> visit = [&](MaterialNodeId id) {
        YBI_ASSERT(id < graph.nodes.size());
        if (reachable[id])
            return;

        reachable[id] = true;
        for (const MaterialInput &input : graph.nodes[id].inputs)
        {
            if (const MaterialNodeId *child = std::get_if<MaterialNodeId>(&input.value))
                visit(*child);
        }
    };

    visit(graph.outputNode);
    return reachable;
}

static MaterialCompilationContext ConvertGraph(
    const MaterialGraph &graph,
    const std::unordered_map<const MaterialNode *, MaterialTextureHandle> &textureNodeToHandle)
{
    MaterialCompilationContext context(textureNodeToHandle);
    YBI_ASSERT(graph.outputNode < graph.nodes.size());

    // Shared default for every unconnected MaterialX normal input. An unused node is
    // excluded from compilationOrder by the traversal from graph.outputNode.
    MaterialNodeBase *defaultNormal =
        context.AddNode<ExpressionNode>(MaterialOpcode::NormalVector3);

    const std::vector<bool> reachable = FindReachableNodes(graph);
    SourceNodeMap sourceNodes(graph.nodes.size(), nullptr);

    // Create every source node before resolving any connections.
    for (MaterialNodeId oldId = 0; oldId < graph.nodes.size(); ++oldId)
    {
        if (!reachable[oldId])
            continue;

        const MaterialNode &node = graph.nodes[oldId];
        MaterialNodeBase *newNode = nullptr;

        switch (node.type)
        {
            case MaterialOpcode::ConstantFloat:
            {
                const MaterialInput *input = FindInput(node, "value");
                const float *value = input ? std::get_if<float>(&input->value) : nullptr;
                YBI_ASSERT(value);
                newNode = context.AddNode<ConstantFloatNode>(*value);
                break;
            }
            case MaterialOpcode::ConstantColor:
            {
                const MaterialInput *input = FindInput(node, "value");
                const Vec3f *value = input ? std::get_if<Vec3f>(&input->value) : nullptr;
                YBI_ASSERT(value);
                newNode = context.AddNode<ConstantColor3Node>(*value);
                break;
            }
            case MaterialOpcode::ConstantVector2:
            {
                const MaterialInput *input = FindInput(node, "value");
                const Vec2f *value = input ? std::get_if<Vec2f>(&input->value) : nullptr;
                YBI_ASSERT(value);
                newNode = context.AddNode<ConstantVector2Node>(*value);
                break;
            }
            case MaterialOpcode::RemapVector3:
            case MaterialOpcode::RemapFloat:
            {
                newNode = context.AddNode<RemapNode>(GetMaterialOpcodeOutputType(node.type));
                break;
            }
            case MaterialOpcode::ArtisticIorIor:
                newNode = context.AddNode<ArtisticIorIorNode>();
                break;
            case MaterialOpcode::ArtisticIorExtinction:
                newNode = context.AddNode<ArtisticIorExtinctionNode>();
                break;
            case MaterialOpcode::NormalMapFloat:
                newNode = context.AddNode<NormalMapNode>();
                break;
            case MaterialOpcode::ViewDirectionVector3:
            case MaterialOpcode::NormalVector3:
            case MaterialOpcode::DotProductVector3:
            case MaterialOpcode::AbsValFloat:
            case MaterialOpcode::AsinFloat:
            case MaterialOpcode::PowerFloat:
            case MaterialOpcode::DivideFloat:
            case MaterialOpcode::DivideVector3:
            case MaterialOpcode::MultiplyFloat:
            case MaterialOpcode::MultiplyVector3:
            case MaterialOpcode::ExpVector3:
            case MaterialOpcode::LnVector3:
            case MaterialOpcode::NormalizeVector3:
            case MaterialOpcode::ConvertFloatVector2:
            case MaterialOpcode::ConvertFloatVector3:
            case MaterialOpcode::ConvertVector3Color3:
            case MaterialOpcode::ConvertColor3Vector3:
            case MaterialOpcode::ExtractVector3:
            case MaterialOpcode::AddFloat:
            case MaterialOpcode::AddVector3:
            case MaterialOpcode::InvertFloat:
            case MaterialOpcode::MixFloat:
            case MaterialOpcode::MixVector2:
            case MaterialOpcode::MixVector3:
            case MaterialOpcode::ContrastFloat:
            case MaterialOpcode::ClampFloat:
                newNode = context.AddNode<ExpressionNode>(node.type);
                break;
            case MaterialOpcode::ImageFloat:
            case MaterialOpcode::ImageVector2:
            case MaterialOpcode::ImageVector3:
            case MaterialOpcode::ImageVector4:
            case MaterialOpcode::TiledImageFloat:
            case MaterialOpcode::TiledImageVector2:
            case MaterialOpcode::TiledImageVector3:
            case MaterialOpcode::TiledImageVector4:
                newNode =
                    context.AddNode<ImageNode>(GetMaterialOpcodeOutputType(node.type), &node);
                break;
            case MaterialOpcode::sheen_bsdf:
                newNode = context.AddNode<SheenBsdfNode>();
                break;
            case MaterialOpcode::dielectric_bsdf:
                newNode = context.AddNode<DielectricBsdfNode>();
                break;
            case MaterialOpcode::oren_nayar_diffuse_bsdf:
                newNode = context.AddNode<OrenNayarBsdfNode>();
                break;
            case MaterialOpcode::conductor_bsdf:
                newNode = context.AddNode<ConductorBsdfNode>();
                break;
            case MaterialOpcode::mix_bsdf:
                newNode = context.AddNode<MixBsdfNode>();
                break;
            case MaterialOpcode::layer_bsdf:
                newNode = context.AddNode<LayerBsdfNode>();
                break;
            default:
                YBI_LOGFATAL("Material compilation does not support opcode %u yet\n",
                             static_cast<uint32_t>(node.type));
                break;
        }

        sourceNodes[oldId] = newNode;
    }

    // Resolve connections and materialize runtime literals as constant nodes.
    for (MaterialNodeId oldId = 0; oldId < graph.nodes.size(); ++oldId)
    {
        if (!reachable[oldId])
            continue;

        const MaterialNode &node = graph.nodes[oldId];
        MaterialNodeBase *newNode = GetMappedNode(sourceNodes, oldId);

        switch (node.type)
        {
            case MaterialOpcode::TiledImageFloat:
            case MaterialOpcode::TiledImageVector2:
            case MaterialOpcode::TiledImageVector3:
            case MaterialOpcode::TiledImageVector4:
            {
                YBI_ASSERT(!FindInput(node, "texcoord") &&
                           "Explicit tiledimage texcoord requires UV derivatives");
                auto *image = static_cast<ImageNode *>(newNode);
                image->uvTiling =
                    AddVector2Input(context, sourceNodes, node, "uvtiling", Vec2f(1.f));
                image->uvOffset =
                    AddVector2Input(context, sourceNodes, node, "uvoffset", Vec2f(0.f));
                image->realWorldImageSize =
                    AddVector2Input(context, sourceNodes, node, "realworldimagesize", Vec2f(1.f));
                image->realWorldTileSize =
                    AddVector2Input(context, sourceNodes, node, "realworldtilesize", Vec2f(1.f));
                image->inputs = {image->uvTiling,
                                 image->uvOffset,
                                 image->realWorldImageSize,
                                 image->realWorldTileSize};
                break;
            }
            case MaterialOpcode::ViewDirectionVector3:
            case MaterialOpcode::NormalVector3:
                AssertStringInput(node, "space", "world");
                break;
            case MaterialOpcode::DotProductVector3:
            case MaterialOpcode::DivideVector3:
            case MaterialOpcode::MultiplyVector3:
            case MaterialOpcode::AddVector3:
            {
                auto *expression = static_cast<ExpressionNode *>(newNode);
                expression->inputs = {
                    AddColor3Input(context, sourceNodes, node, "in1", Vec3f(0.f)),
                    AddColor3Input(context,
                                   sourceNodes,
                                   node,
                                   "in2",
                                   node.type == MaterialOpcode::DivideVector3 ? Vec3f(1.f)
                                                                              : Vec3f(0.f))};
                break;
            }
            case MaterialOpcode::DivideFloat:
            case MaterialOpcode::AddFloat:
            case MaterialOpcode::PowerFloat:
            case MaterialOpcode::MultiplyFloat:
            {
                auto *expression = static_cast<ExpressionNode *>(newNode);
                expression->inputs = {AddFloatInput(context, sourceNodes, node, "in1", 0.f),
                                      AddFloatInput(context,
                                                    sourceNodes,
                                                    node,
                                                    "in2",
                                                    (node.type == MaterialOpcode::DivideFloat ||
                                                     node.type == MaterialOpcode::PowerFloat ||
                                                     node.type == MaterialOpcode::MultiplyFloat)
                                                        ? 1.f
                                                        : 0.f)};
                break;
            }
            case MaterialOpcode::InvertFloat:
            {
                static_cast<ExpressionNode *>(newNode)->inputs = {
                    AddFloatInput(context, sourceNodes, node, "in", 0.f),
                    AddFloatInput(context, sourceNodes, node, "amount", 1.f)};
                break;
            }
            case MaterialOpcode::AbsValFloat:
            case MaterialOpcode::AsinFloat:
            case MaterialOpcode::ConvertFloatVector2:
            case MaterialOpcode::ConvertFloatVector3:
            {
                static_cast<ExpressionNode *>(newNode)->inputs = {
                    AddFloatInput(context, sourceNodes, node, "in", 0.f)};
                break;
            }
            case MaterialOpcode::ExpVector3:
            case MaterialOpcode::LnVector3:
            case MaterialOpcode::NormalizeVector3:
            case MaterialOpcode::ConvertVector3Color3:
            case MaterialOpcode::ConvertColor3Vector3:
            {
                static_cast<ExpressionNode *>(newNode)->inputs = {
                    AddColor3Input(context, sourceNodes, node, "in", Vec3f(0.f))};
                break;
            }
            case MaterialOpcode::ExtractVector3:
            {
                auto *expression = static_cast<ExpressionNode *>(newNode);
                expression->inputs = {
                    AddColor3Input(context, sourceNodes, node, "in", Vec3f(0.f))};
                const MaterialInput *index = FindInput(node, "index");
                if (index)
                {
                    const int *component = std::get_if<int>(&index->value);
                    YBI_ASSERT(component && *component >= 0 && *component < 3);
                    expression->component = *component;
                }
                break;
            }
            case MaterialOpcode::MixFloat:
            {
                static_cast<ExpressionNode *>(newNode)->inputs = {
                    AddFloatInput(context, sourceNodes, node, "bg", 0.f),
                    AddFloatInput(context, sourceNodes, node, "fg", 0.f),
                    AddFloatInput(context, sourceNodes, node, "mix", 0.f)};
                break;
            }
            case MaterialOpcode::MixVector3:
            {
                static_cast<ExpressionNode *>(newNode)->inputs = {
                    AddColor3Input(context, sourceNodes, node, "bg", Vec3f(0.f)),
                    AddColor3Input(context, sourceNodes, node, "fg", Vec3f(0.f)),
                    AddFloatInput(context, sourceNodes, node, "mix", 0.f)};
                break;
            }
            case MaterialOpcode::MixVector2:
            {
                static_cast<ExpressionNode *>(newNode)->inputs = {
                    AddVector2Input(context, sourceNodes, node, "bg", Vec2f(0.f)),
                    AddVector2Input(context, sourceNodes, node, "fg", Vec2f(0.f)),
                    AddFloatInput(context, sourceNodes, node, "mix", 0.f)};
                break;
            }
            case MaterialOpcode::ContrastFloat:
            {
                static_cast<ExpressionNode *>(newNode)->inputs = {
                    AddFloatInput(context, sourceNodes, node, "in", 0.f),
                    AddFloatInput(context, sourceNodes, node, "amount", 1.f),
                    AddFloatInput(context, sourceNodes, node, "pivot", 0.5f)};
                break;
            }
            case MaterialOpcode::ClampFloat:
            {
                static_cast<ExpressionNode *>(newNode)->inputs = {
                    AddFloatInput(context, sourceNodes, node, "in", 0.f),
                    AddFloatInput(context, sourceNodes, node, "low", 0.f),
                    AddFloatInput(context, sourceNodes, node, "high", 1.f)};
                break;
            }
            case MaterialOpcode::RemapFloat:
            {
                auto *remap = static_cast<RemapNode *>(newNode);
                remap->value = AddFloatInput(context, sourceNodes, node, "in", 0.f);
                remap->inLow = AddFloatInput(context, sourceNodes, node, "inlow", 0.f);
                remap->inHigh = AddFloatInput(context, sourceNodes, node, "inhigh", 1.f);
                remap->outLow = AddFloatInput(context, sourceNodes, node, "outlow", 0.f);
                remap->outHigh = AddFloatInput(context, sourceNodes, node, "outhigh", 1.f);
                remap->inputs = {
                    remap->value, remap->inLow, remap->inHigh, remap->outLow, remap->outHigh};
                break;
            }
            case MaterialOpcode::RemapVector3:
            {
                auto *remap = static_cast<RemapNode *>(newNode);
                remap->value = AddColor3Input(context, sourceNodes, node, "in", Vec3f(0.f));
                remap->inLow = AddColor3Input(context, sourceNodes, node, "inlow", Vec3f(0.f));
                remap->inHigh = AddColor3Input(context, sourceNodes, node, "inhigh", Vec3f(1.f));
                remap->outLow = AddColor3Input(context, sourceNodes, node, "outlow", Vec3f(0.f));
                remap->outHigh = AddColor3Input(context, sourceNodes, node, "outhigh", Vec3f(1.f));
                remap->inputs = {
                    remap->value, remap->inLow, remap->inHigh, remap->outLow, remap->outHigh};
                break;
            }
            case MaterialOpcode::ArtisticIorIor:
            {
                auto *artisticIor = static_cast<ArtisticIorIorNode *>(newNode);
                artisticIor->reflectivity = AddColor3Input(
                    context, sourceNodes, node, "reflectivity", Vec3f(0.944f, 0.776f, 0.373f));
                artisticIor->edgeColor = AddColor3Input(
                    context, sourceNodes, node, "edge_color", Vec3f(0.998f, 0.981f, 0.751f));
                artisticIor->inputs = {artisticIor->reflectivity, artisticIor->edgeColor};
                break;
            }
            case MaterialOpcode::ArtisticIorExtinction:
            {
                auto *artisticIor = static_cast<ArtisticIorExtinctionNode *>(newNode);
                artisticIor->reflectivity = AddColor3Input(
                    context, sourceNodes, node, "reflectivity", Vec3f(0.944f, 0.776f, 0.373f));
                artisticIor->edgeColor = AddColor3Input(
                    context, sourceNodes, node, "edge_color", Vec3f(0.998f, 0.981f, 0.751f));
                artisticIor->inputs = {artisticIor->reflectivity, artisticIor->edgeColor};
                break;
            }
            case MaterialOpcode::NormalMapFloat:
            {
                AssertStringInput(node, "space", "tangent");
                auto *normalMap = static_cast<NormalMapNode *>(newNode);
                normalMap->value =
                    AddColor3Input(context, sourceNodes, node, "in", Vec3f(0.5f, 0.5f, 1.f));
                normalMap->scale = AddFloatInput(context, sourceNodes, node, "scale", 1.f);
                normalMap->normal = AddNormalInput(context, sourceNodes, node, defaultNormal);
                normalMap->inputs = {normalMap->value, normalMap->scale, normalMap->normal};
                if (FindInput(node, "tangent"))
                {
                    normalMap->tangent =
                        AddColor3Input(context, sourceNodes, node, "tangent", Vec3f(0.f));
                    normalMap->inputs.push_back(normalMap->tangent);
                }
                if (FindInput(node, "bitangent"))
                {
                    normalMap->bitangent =
                        AddColor3Input(context, sourceNodes, node, "bitangent", Vec3f(0.f));
                    normalMap->inputs.push_back(normalMap->bitangent);
                }
                break;
            }
            case MaterialOpcode::sheen_bsdf:
            {
                if (const MaterialInput *mode = FindInput(node, "mode"))
                {
                    const std::string *modeName = std::get_if<std::string>(&mode->value);
                    YBI_ASSERT(modeName && *modeName == "conty_kulla");
                }

                SheenBsdfNode *sheen = static_cast<SheenBsdfNode *>(newNode);
                sheen->weight = AddFloatInput(context, sourceNodes, node, "weight", 1.0f);
                sheen->color = AddColor3Input(context, sourceNodes, node, "color", Vec3f(1.0f));
                sheen->roughness = AddFloatInput(context, sourceNodes, node, "roughness", 0.3f);
                sheen->inputs = {sheen->weight, sheen->color, sheen->roughness};
                sheen->normal = AddNormalInput(context, sourceNodes, node, defaultNormal);
                sheen->inputs.push_back(sheen->normal);
                break;
            }
            case MaterialOpcode::dielectric_bsdf:
            {
                AssertStringInput(node, "distribution", "ggx");
                AssertStringInput(node, "scatter_mode", "R");
                AssertZeroFloatInput(node, "thinfilm_thickness");

                DielectricBsdfNode *dielectric = static_cast<DielectricBsdfNode *>(newNode);
                dielectric->weight = AddFloatInput(context, sourceNodes, node, "weight", 1.0f);
                dielectric->tint = AddColor3Input(context, sourceNodes, node, "tint", Vec3f(1.0f));
                dielectric->ior = AddFloatInput(context, sourceNodes, node, "ior", 1.5f);
                dielectric->roughness =
                    AddVector2Input(context, sourceNodes, node, "roughness", Vec2f(0.05f));
                dielectric->inputs = {
                    dielectric->weight, dielectric->tint, dielectric->ior, dielectric->roughness};
                dielectric->normal = AddNormalInput(context, sourceNodes, node, defaultNormal);
                dielectric->inputs.push_back(dielectric->normal);
                break;
            }
            case MaterialOpcode::oren_nayar_diffuse_bsdf:
            {
                AssertFalseBoolInput(node, "energy_compensation");

                OrenNayarBsdfNode *diffuse = static_cast<OrenNayarBsdfNode *>(newNode);
                diffuse->weight = AddFloatInput(context, sourceNodes, node, "weight", 1.0f);
                diffuse->color = AddColor3Input(context, sourceNodes, node, "color", Vec3f(0.18f));
                diffuse->roughness = AddFloatInput(context, sourceNodes, node, "roughness", 0.0f);
                diffuse->inputs = {diffuse->weight, diffuse->color, diffuse->roughness};
                diffuse->normal = AddNormalInput(context, sourceNodes, node, defaultNormal);
                diffuse->inputs.push_back(diffuse->normal);
                break;
            }
            case MaterialOpcode::conductor_bsdf:
            {
                AssertStringInput(node, "distribution", "ggx");
                AssertZeroFloatInput(node, "thinfilm_thickness");

                auto *conductor = static_cast<ConductorBsdfNode *>(newNode);
                conductor->weight = AddFloatInput(context, sourceNodes, node, "weight", 1.f);
                conductor->ior =
                    AddColor3Input(context, sourceNodes, node, "ior", Vec3f(0.18f, 0.42f, 1.37f));
                conductor->extinction = AddColor3Input(
                    context, sourceNodes, node, "extinction", Vec3f(3.42f, 2.35f, 1.77f));
                conductor->roughness =
                    AddVector2Input(context, sourceNodes, node, "roughness", Vec2f(0.f));
                conductor->inputs = {conductor->weight,
                                     conductor->ior,
                                     conductor->extinction,
                                     conductor->roughness};
                conductor->normal = AddNormalInput(context, sourceNodes, node, defaultNormal);
                conductor->inputs.push_back(conductor->normal);
                break;
            }
            case MaterialOpcode::mix_bsdf:
            {
                auto *mix = static_cast<MixBsdfNode *>(newNode);
                mix->bg = GetNodeInput(sourceNodes, node, "bg");
                mix->fg = GetNodeInput(sourceNodes, node, "fg");
                mix->mix = AddFloatInput(context, sourceNodes, node, "mix", 0.f);
                mix->inputs = {mix->bg, mix->fg, mix->mix};
                break;
            }
            case MaterialOpcode::layer_bsdf:
            {
                LayerBsdfNode *layer = static_cast<LayerBsdfNode *>(newNode);
                layer->top = GetNodeInput(sourceNodes, node, "top");
                layer->base = GetNodeInput(sourceNodes, node, "base");
                layer->inputs = {layer->top, layer->base};
                break;
            }
            default:
                break;
        }
    }

    MaterialNodeBase *outputNode = GetMappedNode(sourceNodes, graph.outputNode);
    context.TopologicalSort(outputNode);

    return context;
}

} // namespace

std::string MaterialCompiler::Compile(const MaterialGraph &graph)
{
    static const std::unordered_map<const MaterialNode *, MaterialTextureHandle> emptyTextureMap;
    return Compile(graph, emptyTextureMap);
}

std::string MaterialCompiler::Compile(
    const MaterialGraph &graph,
    const std::unordered_map<const MaterialNode *, MaterialTextureHandle> &textureNodeToHandle)
{
    MaterialCompilationContext context = ConvertGraph(graph, textureNodeToHandle);
    const uint32_t materialIndex = numMaterials++;
    const std::string materialPrefix = "m" + std::to_string(materialIndex);

    bool hasSheen = false;
    bool hasDielectric = false;
    bool hasOrenNayar = false;
    bool hasConductor = false;
    bool hasImage = false;
    for (MaterialNodeBase *node : context.compilationOrder)
    {
        hasSheen |= node->GetNodeType() == MaterialOpcode::sheen_bsdf;
        hasDielectric |= node->GetNodeType() == MaterialOpcode::dielectric_bsdf;
        hasOrenNayar |= node->GetNodeType() == MaterialOpcode::oren_nayar_diffuse_bsdf;
        hasConductor |= node->GetNodeType() == MaterialOpcode::conductor_bsdf;
        hasImage |= IsTextureOpcode(node->GetNodeType());
    }

    std::string result;
    if (hasSheen)
        result += "#include \"kernel/bsdf/materialx/sheen.h\"\n";
    if (hasDielectric)
        result += "#include \"kernel/bsdf/materialx/dielectric_bsdf.h\"\n";
    if (hasOrenNayar)
        result += "#include \"kernel/bsdf/materialx/oren_nayar_bsdf.h\"\n";
    if (hasConductor)
        result += "#include \"kernel/bsdf/materialx/conductor_bsdf.h\"\n";
    if (hasImage)
        result += "#include \"kernel/shading/image.h\"\n";
    result += "#include \"kernel/bsdf/materialx/bsdf_frame.h\"\n";
    result += "\nnamespace ybi\n{\n\n";

    std::unordered_set<const MaterialNodeBase *> closureInputs;
    for (const MaterialNodeBase *node : context.compilationOrder)
    {
        if (!node->IsBsdf())
            continue;

        for (const MaterialNodeBase *input : node->inputs)
        {
            if (!input->IsBsdf() && !IsConstantOpcode(input->GetNodeType()))
                closureInputs.insert(input);
        }
    }

    size_t closureOffset = 0;
    for (MaterialNodeBase *node : context.compilationOrder)
    {
        result += node->Compile(context, materialPrefix);

        if (closureInputs.count(node) != 0)
        {
            result += node->CompileClosureAccessors(materialPrefix, closureOffset);
            closureOffset += GetMaterialOpcodeOutputRegisterCount(node->GetNodeType());
        }

        if (node->IsBsdf())
        {
            result +=
                static_cast<BsdfNodeBase *>(node)->CompileInitialization(context, materialPrefix);
        }
    }

    YBI_ASSERT(closureOffset <= MAX_CLOSURE_SIZE);

    result += "} // namespace ybi\n";
    generatedMaterialSources.push_back(result);
    return result;
}

MaterialCompilationResult MaterialCompiler::Finalize(Device &device) const
{
    MaterialCompilationResult result;
    result.source = "#include \"kernel/launch/kernel_params.h\"\n"
                    "#include \"kernel/bsdf/closure_types.h\"\n\n";
#ifdef YBI_MATERIAL_DEBUG
    result.source += "extern \"C\" __device__ int printf(const char *, ...);\n\n";
#endif

    for (const std::string &source : generatedMaterialSources)
    {
        result.source += source;
        result.source += "\n";
    }

    std::vector<std::string> entryNames;
    entryNames.reserve(generatedMaterialSources.size() * 4);

    const std::string callableSignature = device.GetCallableFunctionSignature();
    const std::string directCallablePrefix = device.GetDirectCallablePrefix();

    result.source += "namespace ybi\n{\n\n";
    for (uint32_t materialIndex = 0;
         materialIndex < static_cast<uint32_t>(generatedMaterialSources.size());
         ++materialIndex)
    {
        const std::string materialName = "m" + std::to_string(materialIndex);
        const std::string rootName = materialName + "_n0";

        const std::string initializeName = INITIALIZE_PREFIX + materialName;
        const std::string sampleName = SAMPLE_PREFIX + materialName;
        const std::string evaluateName = EVALUATE_PREFIX + materialName;
        const std::string pdfName = PDF_PREFIX + materialName;

        const auto compileCallable = [&](const char *stubPath,
                                         const std::string &name,
                                         const std::string &implementation) {
            return CompileWithMacros(stubPath,
                                     {{"YBI_MATERIAL_CALLABLE_SIGNATURE", callableSignature},
                                      {"YBI_MATERIAL_CALLABLE_NAME", directCallablePrefix + name},
                                      {"YBI_MATERIAL_CALLABLE_IMPL", implementation}});
        };

        result.source += compileCallable("kernel/bsdf/materialx/stubs/initialize_callable_impl.h",
                                         initializeName,
                                         INITIALIZE_PREFIX + rootName);
        result.source += compileCallable("kernel/bsdf/materialx/stubs/sample_callable_impl.h",
                                         sampleName,
                                         SAMPLE_PREFIX + rootName);
        result.source += compileCallable("kernel/bsdf/materialx/stubs/evaluate_callable_impl.h",
                                         evaluateName,
                                         EVALUATE_PREFIX + rootName);
        result.source += compileCallable(
            "kernel/bsdf/materialx/stubs/pdf_callable_impl.h", pdfName, PDF_PREFIX + rootName);

        entryNames.push_back(directCallablePrefix + initializeName);
        entryNames.push_back(directCallablePrefix + sampleName);
        entryNames.push_back(directCallablePrefix + evaluateName);
        entryNames.push_back(directCallablePrefix + pdfName);
    }
    result.source += "} // namespace ybi\n";

    if (entryNames.empty())
    {
        result.succeeded = true;
        return result;
    }

    DeviceModule module = 0;
    if (!device.CompileKernel(result.source, module))
        return result;

    std::vector<uint32_t> callableIndices;
    if (!device.RegisterDirectCallables(module, entryNames, callableIndices))
        return result;

    YBI_ASSERT(callableIndices.size() == entryNames.size());
    result.materials.resize(generatedMaterialSources.size());
    for (uint32_t materialIndex = 0;
         materialIndex < static_cast<uint32_t>(result.materials.size());
         ++materialIndex)
    {
        DeviceMaterial3 &material = result.materials[materialIndex];
        material.initializeCallableIndex = callableIndices[4 * materialIndex + 0];
        material.sampleCallableIndex = callableIndices[4 * materialIndex + 1];
        material.evaluateCallableIndex = callableIndices[4 * materialIndex + 2];
        material.pdfCallableIndex = callableIndices[4 * materialIndex + 3];
    }

    result.succeeded = true;
    return result;
}

} // namespace ybi
