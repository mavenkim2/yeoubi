#include "scene/material_compilation.h"
#include "core/assert.h"
#include "device/device.h"
#include "render/bsdf/materialx/bsdf_types.h"
#include "render/bsdf/materialx/dielectric_bsdf.h"
#include <functional>
#include <initializer_list>
#include <iomanip>
#include <limits>
#include <memory>
#include <sstream>
#include <unordered_map>
#include <utility>
#include <vector>

namespace ybi
{
namespace
{

static const std::string GET_CLOSURE_PREFIX = "GetClosure_";

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
        result += "#undef " + std::string(definition.first) + "\n";

    return result;
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
                "render/bsdf/materialx/stubs/closure_impl.h",
                {
                    {"YBI_MATERIAL_VALUE_TYPE", returnType},
                    {"YBI_MATERIAL_SELF_GET_CLOSURE", GET_CLOSURE_PREFIX + nodeName},
                    {"YBI_MATERIAL_SELF_SET_CLOSURE", "SetClosure_" + nodeName},
                    {"YBI_MATERIAL_VALUE", "Value_" + nodeName},
                    {"YBI_MATERIAL_CLOSURE_OFFSET", std::to_string(closureOffset)},
                });
        }
        return "";
    }

    uint32_t index = kInvalidMaterialNodeId;
    std::vector<MaterialNodeBase *> inputs;
};

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
        std::string source = "YBI_GPU void Initialize_" + FullName(materialPrefix) +
                             "(const KernelParams &params, ShadingInfo &info)\n{\n    ";

        for (MaterialNodeBase *node : inputs)
        {
            if (node->IsBsdf())
            {
                source += "Initialize_" + node->FullName(materialPrefix) + "(params, info);\n";
            }
            else if (!node->IsBsdf() && !IsConstantOpcode(node->GetNodeType()))
            {
                std::string nodeName = node->FullName(materialPrefix);
                const char *returnType = GetMaterialOpcodeOutputString(node->GetNodeType());

                source += "SetClosure_" + nodeName + "(params, info);\n";
            }
        }
        source += "}\n";
        return source;
    }
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
        return "YBI_GPU float " + GET_CLOSURE_PREFIX + FullName(materialPrefix) +
               "(const ShadingInfo &)\n{\n    return " + FloatLiteral(value) + ";\n}\n";
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
        return "YBI_GPU Vec3f " + GET_CLOSURE_PREFIX + FullName(materialPrefix) +
               "(const ShadingInfo &)\n{\n    return Vec3f(" + FloatLiteral(value.x) + ", " +
               FloatLiteral(value.y) + ", " + FloatLiteral(value.z) + ");\n}\n";
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
        return "YBI_GPU Vec2f " + GET_CLOSURE_PREFIX + FullName(materialPrefix) +
               "(const ShadingInfo &)\n{\n    return Vec2f(" + FloatLiteral(value.x) + ", " +
               FloatLiteral(value.y) + ");\n}\n";
    }

    MaterialOpcode GetNodeType() const override
    {
        return MaterialOpcode::ConstantVector2;
    }

    Vec2f value;
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

        return CompileWithMacros("render/bsdf/materialx/stubs/image_impl.h",
                                 {{"YBI_MATERIAL_RETURN", returnType},
                                  {"YBI_MATERIAL_SELF", "Value_" + FullName(materialPrefix)},
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

        return CompileWithMacros("render/bsdf/materialx/stubs/sheen_impl.h",
                                 {{"YBI_MATERIAL_SELF_SAMPLE", "Sample_" + self},
                                  {"YBI_MATERIAL_SELF_EVALUATE", "Evaluate_" + self},
                                  {"YBI_MATERIAL_SELF_PDF", "Pdf_" + self},
                                  {"YBI_MATERIAL_SELF_ALBEDO", "Albedo_" + self},
                                  {"YBI_MATERIAL_SELF_TRANSMITTANCE", "Transmittance_" + self},
                                  {"YBI_MATERIAL_WEIGHT", GET_CLOSURE_PREFIX + weightName},
                                  {"YBI_MATERIAL_COLOR", GET_CLOSURE_PREFIX + colorName},
                                  {"YBI_MATERIAL_ROUGHNESS", GET_CLOSURE_PREFIX + roughnessName}});
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

        return CompileWithMacros("render/bsdf/materialx/stubs/dielectric_impl.h",
                                 {{"YBI_MATERIAL_SELF_SAMPLE", "Sample_" + self},
                                  {"YBI_MATERIAL_SELF_EVALUATE", "Evaluate_" + self},
                                  {"YBI_MATERIAL_SELF_PDF", "Pdf_" + self},
                                  {"YBI_MATERIAL_SELF_ALBEDO", "Albedo_" + self},
                                  {"YBI_MATERIAL_SELF_TRANSMITTANCE", "Transmittance_" + self},
                                  {"YBI_MATERIAL_WEIGHT", GET_CLOSURE_PREFIX + weightName},
                                  {"YBI_MATERIAL_TINT", GET_CLOSURE_PREFIX + tintName},
                                  {"YBI_MATERIAL_IOR", GET_CLOSURE_PREFIX + iorName},
                                  {"YBI_MATERIAL_ROUGHNESS", GET_CLOSURE_PREFIX + roughnessName}});
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

        return CompileWithMacros("render/bsdf/materialx/stubs/oren_nayar_impl.h",
                                 {{"YBI_MATERIAL_SELF_SAMPLE", "Sample_" + self},
                                  {"YBI_MATERIAL_SELF_EVALUATE", "Evaluate_" + self},
                                  {"YBI_MATERIAL_SELF_PDF", "Pdf_" + self},
                                  {"YBI_MATERIAL_SELF_ALBEDO", "Albedo_" + self},
                                  {"YBI_MATERIAL_SELF_TRANSMITTANCE", "Transmittance_" + self},
                                  {"YBI_MATERIAL_WEIGHT", GET_CLOSURE_PREFIX + weightName},
                                  {"YBI_MATERIAL_COLOR", GET_CLOSURE_PREFIX + colorName},
                                  {"YBI_MATERIAL_ROUGHNESS", GET_CLOSURE_PREFIX + roughnessName}});
    }

    MaterialOpcode GetNodeType() const override
    {
        return MaterialOpcode::oren_nayar_diffuse_bsdf;
    }

    MaterialNodeBase *weight = nullptr;
    MaterialNodeBase *color = nullptr;
    MaterialNodeBase *roughness = nullptr;
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
            "render/bsdf/materialx/stubs/layer_impl.h",
            {{"YBI_MATERIAL_SELF_SAMPLE", "Sample_" + self},
             {"YBI_MATERIAL_SELF_EVALUATE", "Evaluate_" + self},
             {"YBI_MATERIAL_SELF_PDF", "Pdf_" + self},
             {"YBI_MATERIAL_SELF_ALBEDO", "Albedo_" + self},
             {"YBI_MATERIAL_SELF_TRANSMITTANCE", "Transmittance_" + self},
             {"YBI_MATERIAL_SAMPLE_TOP", "Sample_" + topName},
             {"YBI_MATERIAL_EVALUATE_TOP", "Evaluate_" + topName},
             {"YBI_MATERIAL_PDF_TOP", "Pdf_" + topName},
             {"YBI_MATERIAL_ALBEDO_TOP", "Albedo_" + topName},
             {"YBI_MATERIAL_TRANSMITTANCE_TOP", "Transmittance_" + topName},
             {"YBI_MATERIAL_SAMPLE_BOTTOM", "Sample_" + baseName},
             {"YBI_MATERIAL_EVALUATE_BOTTOM", "Evaluate_" + baseName},
             {"YBI_MATERIAL_PDF_BOTTOM", "Pdf_" + baseName},
             {"YBI_MATERIAL_ALBEDO_BOTTOM", "Albedo_" + baseName},
             {"YBI_MATERIAL_TRANSMITTANCE_BOTTOM", "Transmittance_" + baseName}});
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
            case MaterialOpcode::ImageFloat:
            case MaterialOpcode::ImageVector2:
            case MaterialOpcode::ImageVector3:
            case MaterialOpcode::ImageVector4:
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
    bool hasImage = false;
    for (MaterialNodeBase *node : context.compilationOrder)
    {
        hasSheen |= node->GetNodeType() == MaterialOpcode::sheen_bsdf;
        hasDielectric |= node->GetNodeType() == MaterialOpcode::dielectric_bsdf;
        hasOrenNayar |= node->GetNodeType() == MaterialOpcode::oren_nayar_diffuse_bsdf;
        hasImage |= IsTextureOpcode(node->GetNodeType());
    }

    std::string result;
    if (hasSheen)
        result += "#include \"render/bsdf/materialx/sheen.h\"\n";
    if (hasDielectric)
        result += "#include \"render/bsdf/materialx/dielectric_bsdf.h\"\n";
    if (hasOrenNayar)
        result += "#include \"render/bsdf/materialx/oren_nayar_bsdf.h\"\n";
    if (hasImage)
        result += "#include \"render/shaders/image.h\"\n";
    result += "\nnamespace ybi\n{\n\n";

    size_t closureOffset = 0;
    for (MaterialNodeBase *node : context.compilationOrder)
    {
        result += node->Compile(context, materialPrefix);

        // If the node is used to poopulate the closure, compile those accessors
        if (!node->IsBsdf() && !IsConstantOpcode(node->GetNodeType()))
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
    result.source = "#include \"device/kernel_params.h\"\n"
                    "#include \"render/shaders/shading.h\"\n\n";

    for (const std::string &source : generatedMaterialSources)
    {
        result.source += source;
        result.source += "\n";
    }

    std::vector<std::string> entryNames;
    entryNames.reserve(generatedMaterialSources.size() * 3);

    result.source += "namespace ybi\n{\n\n";
    for (uint32_t materialIndex = 0;
         materialIndex < static_cast<uint32_t>(generatedMaterialSources.size());
         ++materialIndex)
    {
        const std::string materialName = "m" + std::to_string(materialIndex);
        const std::string rootName = materialName + "_n0";

        const std::string initializeName = "Initialize_" + materialName;
        const std::string sampleName = "Sample_" + materialName;
        const std::string evaluateName = "Evaluate_" + materialName;
        const std::string pdfName = "Pdf_" + materialName;

        result.source +=
            device.GetCallableFunctionSignature(
                "void", initializeName, "const KernelParams *params, ShadingInfo *info") +
            "\n{\n    " + "Initialize_" + rootName + "(*params, *info);\n}\n\n";
        result.source += device.GetCallableFunctionSignature(
                             "bool",
                             sampleName,
                             "const KernelParams *params, const ShadingInfo *info, "
                             "Vec3f viewDirection, Vec3f *lightDirection, Vec3f u") +
                         "\n{\n    return Sample_" + rootName +
                         "(*params, *info, viewDirection, *lightDirection, u);\n}\n\n";

        result.source += device.GetCallableFunctionSignature(
                             "Vec3f",
                             evaluateName,
                             "const KernelParams *params, const ShadingInfo *info, "
                             "Vec3f viewDirection, Vec3f lightDirection") +
                         "\n{\n    return Evaluate_" + rootName +
                         "(*params, *info, viewDirection, lightDirection);\n}\n\n";

        result.source += device.GetCallableFunctionSignature(
                             "float",
                             pdfName,
                             "const KernelParams *params, const ShadingInfo *info, "
                             "Vec3f viewDirection, Vec3f lightDirection") +
                         "\n{\n    return Pdf_" + rootName +
                         "(*params, *info, viewDirection, lightDirection);\n}\n\n";

        entryNames.push_back("__direct_callable__" + initializeName);
        entryNames.push_back("__direct_callable__" + sampleName);
        entryNames.push_back("__direct_callable__" + evaluateName);
        entryNames.push_back("__direct_callable__" + pdfName);
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
