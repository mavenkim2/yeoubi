#include "scene/shader_compilation.h"
#include "core/assert.h"
#include "core/bit.h"
#include "core/bit_vector.h"
#include "core/thread_context.h"
#include "render/material_opcodes.h"
#include "render/material_types.h"
#include "render/virtual_texture_manager.h"
#include "scene/device_scene.h"
#include "scene/material.h"
#include <unordered_map>
#include <vector>

namespace ybi
{

struct NodeInfo
{
    std::vector<MaterialNodeId> outNodes;
    int in = 0;
    int out = 0;

    // Used when compiling
    int inputSlotOffset;
    int numInputSlots;

    int outputSlotOffset;
    int numOutputSlots;
};

struct ShaderCompilationContext
{
    const MaterialGraph &graph;
    std::vector<uint32_t> &nodes;
    const std::unordered_map<const MaterialNode *, MaterialTextureHandle> &textureNodeToHandle;
    const std::unordered_map<std::string, AttributeId> &attributeIdMap;

    void AddFloat(float v)
    {
        nodes.push_back(FloatAsUint(v));
    }
    void AddUint(uint32_t v)
    {
        nodes.push_back(v);
    }
    void AddTextureHandle(MaterialTextureHandle handle)
    {
        nodes.push_back(handle.udimBase);
        nodes.push_back(handle.udimTableOffset);
    }
    AttributeId GetAttributeId(const std::string &name)
    {
        auto find = attributeIdMap.find(name);
        YBI_ASSERT(find != attributeIdMap.end());
        return find->second;
    }
};

template <typename Func>
static void IterateMaterialNodeInputs(const MaterialNode &node, Func &&func)
{
    for (const MaterialInput &input : node.inputs)
    {
        std::visit(
            [&](const auto &value) {
                using T = std::decay_t<decltype(value)>;
                if constexpr (std::is_same_v<T, MaterialNodeId>)
                {
                    func(value);
                }
            },
            input.value);
    }
}

template <typename T>
static const T *GetConstantValue(const MaterialNode &node)
{
    YBI_ERROR(node.inputs.size() == 1 && node.inputs[0].name == "value",
              "%s requires exactly one input named value\n",
              ConvertMaterialOpcodeToString(node.type));
    if (node.inputs.size() != 1 || node.inputs[0].name != "value")
    {
        return nullptr;
    }

    const T *value = std::get_if<T>(&node.inputs[0].value);
    YBI_ERROR(value != nullptr,
              "%s has an input with an unexpected type\n",
              ConvertMaterialOpcodeToString(node.type));
    if (value == nullptr)
    {
        return nullptr;
    }

    return value;
}

template <typename Func>
static void IterateNamedInput(const MaterialNode &node, const std::string &name, Func &&func)
{
    for (const MaterialInput &input : node.inputs)
    {
        if (input.name == name)
        {
            std::visit(
                [&](const auto &value) {
                    using T = std::decay_t<decltype(value)>;
                    func(value);
                },
                input.value);
            break;
        }
    }
}

static void CompileConstantFloat(const MaterialNode &node,
                                 const std::vector<NodeInfo> &nodeInfos,
                                 MaterialNodeId nodeId,
                                 ShaderCompilationContext &context)
{
    const float *value = GetConstantValue<float>(node);
    if (value == nullptr)
        return;

    context.nodes.push_back(static_cast<uint32_t>(node.type));
    context.AddFloat(*value);
    context.AddUint(nodeInfos[nodeId].outputSlotOffset);
}

static void CompileConstantBool(const MaterialNode &node,
                                const std::vector<NodeInfo> &,
                                MaterialNodeId,
                                ShaderCompilationContext &context)
{
    YBI_ASSERT(0);
    const bool *value = GetConstantValue<bool>(node);
    if (value == nullptr)
        return;

    context.nodes.push_back(static_cast<uint32_t>(node.type));
    context.nodes.push_back(*value ? 1u : 0u);
}

static void CompileConstantColor(const MaterialNode &node,
                                 const std::vector<NodeInfo> &nodeInfos,
                                 MaterialNodeId nodeId,
                                 ShaderCompilationContext &context)
{
    const Vec3f *value = GetConstantValue<Vec3f>(node);
    if (value == nullptr)
        return;

    context.nodes.push_back(static_cast<uint32_t>(node.type));
    context.AddFloat(value->x);
    context.AddFloat(value->y);
    context.AddFloat(value->z);
    context.AddUint(nodeInfos[nodeId].outputSlotOffset);
}

// -----------------------------------------------------------------------------
// Image node compilation
// -----------------------------------------------------------------------------

static MaterialTextureHandle GetTextureHandle(const MaterialNode &node,
                                              ShaderCompilationContext &context)
{
    const auto found = context.textureNodeToHandle.find(&node);
    YBI_ERROR(found != context.textureNodeToHandle.end(),
              "Image node does not have a material texture handle\n");

    MaterialTextureHandle handle = found->second;
    return handle;
}

// TODO: add address mode
static void CompileImageVector3(const MaterialNode &node,
                                const std::vector<NodeInfo> &nodeInfos,
                                MaterialNodeId nodeId,
                                ShaderCompilationContext &context)
{
    context.nodes.push_back(static_cast<uint32_t>(node.type));
    MaterialTextureHandle handle = GetTextureHandle(node, context);
    context.AddTextureHandle(handle);
    context.AddUint(nodeInfos[nodeId].outputSlotOffset);
}

static void CompileImageFloat(const MaterialNode &node,
                              const std::vector<NodeInfo> &nodeInfos,
                              MaterialNodeId nodeId,
                              ShaderCompilationContext &context)
{
    context.nodes.push_back(static_cast<uint32_t>(node.type));
    MaterialTextureHandle handle = GetTextureHandle(node, context);
    context.AddTextureHandle(handle);
    context.AddUint(nodeInfos[nodeId].outputSlotOffset);
}

static void CompileAttributeValueVector2(const MaterialNode &node,
                                         const std::vector<NodeInfo> &nodeInfos,
                                         MaterialNodeId nodeId,
                                         ShaderCompilationContext &context)
{
    context.nodes.push_back(static_cast<uint32_t>(node.type));
    bool found = false;
    for (const MaterialInput &input : node.inputs)
    {
        if (input.name == "geomprop")
        {
            const std::string *value = std::get_if<std::string>(&input.value);
            YBI_ASSERT(value);
            AttributeId id = context.GetAttributeId(*value);
            context.AddUint(id);
            context.AddUint(nodeInfos[nodeId].outputSlotOffset);
            found = true;
            break;
        }
    }
    YBI_ASSERT(found);
}

static void CompileHeightToNormalVector3(const MaterialNode &node,
                                         const std::vector<NodeInfo> &,
                                         MaterialNodeId,
                                         ShaderCompilationContext &context)
{
    YBI_ASSERT(0);
}

static void CompileInvertFloat(const MaterialNode &node,
                               const std::vector<NodeInfo> &,
                               MaterialNodeId,
                               ShaderCompilationContext &context)
{
    YBI_ASSERT(0);
}

static void CompileExtractVector3(const MaterialNode &node,
                                  const std::vector<NodeInfo> &nodeInfos,
                                  MaterialNodeId nodeId,
                                  ShaderCompilationContext &context)
{
    const MaterialNodeId *inputNodeId = nullptr;
    int index = 0;

    for (const MaterialInput &input : node.inputs)
    {
        if (input.name == "in")
        {
            inputNodeId = std::get_if<MaterialNodeId>(&input.value);
            YBI_ASSERT(inputNodeId);
        }
        else if (input.name == "index")
        {
            const int *value = std::get_if<int>(&input.value);
            YBI_ASSERT(value);
            index = *value;
        }
    }

    YBI_ASSERT(inputNodeId);
    YBI_ASSERT(*inputNodeId < nodeInfos.size());
    YBI_ASSERT(index >= 0 && index < 3);

    context.AddUint(static_cast<uint32_t>(node.type));
    context.AddUint(nodeInfos[*inputNodeId].outputSlotOffset);
    context.AddUint(static_cast<uint32_t>(index));
    context.AddUint(nodeInfos[nodeId].outputSlotOffset);
}

static void CompileCombine3Vector3(const MaterialNode &node,
                                   const std::vector<NodeInfo> &,
                                   MaterialNodeId,
                                   ShaderCompilationContext &context)
{
    YBI_ASSERT(0);
}

static void CompileNormalMapVector2(const MaterialNode &node,
                                    const std::vector<NodeInfo> &,
                                    MaterialNodeId,
                                    ShaderCompilationContext &context)
{
    YBI_ASSERT(0);
}

static void CompileColorCorrectColor3(const MaterialNode &node,
                                      const std::vector<NodeInfo> &nodeInfos,
                                      MaterialNodeId nodeId,
                                      ShaderCompilationContext &context)
{
    context.nodes.push_back(static_cast<uint32_t>(node.type));
    bool found = false;
    for (const MaterialInput &input : node.inputs)
    {
        if (input.name == "in")
        {
            const MaterialNodeId *inputNodeId = std::get_if<MaterialNodeId>(&input.value);
            if (inputNodeId)
            {
                found = true;
                context.AddUint(nodeInfos[*inputNodeId].outputSlotOffset);
            }
            else
            {
                YBI_ASSERT(0);
            }
        }
    }
    YBI_ASSERT(found);

    found = false;
    for (const MaterialInput &input : node.inputs)
    {
        if (input.name == "gamma")
        {
            const float *gammaVal = std::get_if<float>(&input.value);
            if (gammaVal)
            {
                found = true;
                context.AddUint(0);
                context.AddFloat(*gammaVal);
            }
            else
            {
                YBI_ASSERT(0);
            }
        }
    }
    YBI_ASSERT(found);

    const NodeInfo &nodeInfo = nodeInfos[nodeId];
    context.AddUint(nodeInfo.outputSlotOffset);
}

static void CompileRemapColor3(const MaterialNode &node,
                               const std::vector<NodeInfo> &,
                               MaterialNodeId,
                               ShaderCompilationContext &context)
{
    YBI_ASSERT(0);
}

static void CompileRamp4Float(const MaterialNode &node,
                              const std::vector<NodeInfo> &,
                              MaterialNodeId,
                              ShaderCompilationContext &context)
{
    YBI_ASSERT(0);
}

static void CompileSeparate3Color3(const MaterialNode &node,
                                   const std::vector<NodeInfo> &,
                                   MaterialNodeId,
                                   ShaderCompilationContext &context)
{
    YBI_ASSERT(0);
}

DeviceMaterialExpression CompileShader(
    const MaterialGraph &graph,
    MaterialNodeId outputNode,
    const std::unordered_map<const MaterialNode *, MaterialTextureHandle> &textureNodeToHandle,
    const std::unordered_map<std::string, AttributeId> &attributeIdMap,
    std::vector<uint32_t> &nodes)
{
    ShaderCompilationContext context{graph, nodes, textureNodeToHandle, attributeIdMap};
    const uint32_t bytecodeOffset = static_cast<uint32_t>(nodes.size());

    // Count number of incoming and outgoing edges
    YBI_ASSERT(outputNode < graph.nodes.size());

    std::vector<NodeInfo> nodeInfo(graph.nodes.size());
    std::vector<bool> discovered(graph.nodes.size(), false);
    std::vector<MaterialNodeId> discoveryOrder;
    std::vector<MaterialNodeId> stack;

    stack.push_back(outputNode);
    discovered[outputNode] = true;
    discoveryOrder.push_back(outputNode);

    while (!stack.empty())
    {
        MaterialNodeId node = stack.back();
        stack.pop_back();
        IterateMaterialNodeInputs(graph.nodes[node], [&](MaterialNodeId parent) {
            YBI_ASSERT(parent < graph.nodes.size());

            nodeInfo[node].in++;
            nodeInfo[parent].outNodes.push_back(node);
            nodeInfo[parent].out++;

            if (!discovered[parent])
            {
                discovered[parent] = true;
                discoveryOrder.push_back(parent);
                stack.push_back(parent);
            }
        });
    }

    // Topological sort the nodes
    std::vector<MaterialNodeId> queue;
    for (MaterialNodeId node : discoveryOrder)
    {
        if (nodeInfo[node].in == 0)
        {
            queue.push_back(node);
        }
    }

    std::vector<MaterialNodeId> outputList;
    size_t index = 0;
    while (index != queue.size())
    {
        MaterialNodeId node = queue[index++];
        outputList.push_back(node);
        for (MaterialNodeId neighbor : nodeInfo[node].outNodes)
        {
            YBI_ASSERT(nodeInfo[neighbor].in);
            if (--nodeInfo[neighbor].in == 0)
            {
                queue.push_back(neighbor);
            }
        }
    }

    // Allocate register slots
    ScratchArena scratch;
    BitVector slotAllocated(scratch.arena, SHADER_NUM_REGISTERS);

    // std::vector<uint8_t> compiledBytes;
    // std::vector<CompileNode> compiledNodes;

    for (MaterialNodeId nodeId : outputList)
    {
        const MaterialNode &node = graph.nodes[nodeId];
        size_t numWords = GetMaterialOpcodeOutputRegisterCount(node.type);

#if 0
        // See if node computes a constant
        bool constant = true;
        for (const MaterialInput &input : node.inputs)
        {
            std::visit(
                [&](const auto &value) {
                    using T = std::decay_t<decltype(value)>;

                    if constexpr (std::is_same_v<T, MaterialNodeId>)
                    {
                        constant = false;
                    }
                },
                input.value);
            if (!constant)
                break;
        }

        if (constant)
        {
            printf("%s\n", ConvertMaterialOpcodeToString(node.type));
        }
#endif

        // Allocate output slot
        int freeSlot = 0;
        for (int i = 0; i < SHADER_NUM_REGISTERS; i++)
        {
            if (slotAllocated.GetBit(i))
            {
                freeSlot = i + 1;
            }
            if (i - freeSlot + 1 == numWords)
            {
                for (int i = freeSlot; i < freeSlot + numWords; i++)
                {
                    YBI_ASSERT(!slotAllocated.GetBit(i));
                    slotAllocated.SetBit(i);
                }
                break;
            }
        }

        YBI_ERROR(freeSlot != SHADER_NUM_REGISTERS,
                  "Shader %s ran out of registers during compilation\n",
                  ConvertMaterialOpcodeToString(graph.nodes[outputNode].type));

        // Free slots from any parents whose output is no longer used
        IterateMaterialNodeInputs(node, [&](MaterialNodeId parent) {
            NodeInfo &info = nodeInfo[parent];
            if (--info.out == 0)
            {
                for (int i = info.outputSlotOffset;
                     i < info.outputSlotOffset + info.numOutputSlots;
                     i++)
                {
                    YBI_ASSERT(slotAllocated.GetBit(i));
                    slotAllocated.UnsetBit(i);
                }
            }
        });

        NodeInfo &info = nodeInfo[nodeId];
        info.outputSlotOffset = freeSlot;
        info.numOutputSlots = numWords;
    }

    for (MaterialNodeId nodeId : outputList)
    {
        const MaterialNode &node = graph.nodes[nodeId];
        switch (node.type)
        {
#define CASE(name, type) \
    case MaterialOpcode::name: \
    { \
        Compile##name(node, nodeInfo, nodeId, context); \
        break; \
    }
            RENDERER_NODE_DEFINITIONS(CASE)
#undef CASE
            default:
                YBI_LOGFATAL("Shader compilation for %s is not implemented\n",
                             ConvertMaterialOpcodeToString(node.type));
                break;
        }
    }

    context.AddUint(static_cast<uint32_t>(MaterialOpcode::Count));

    DeviceMaterialExpression result;
    result.bytecodeOffset = bytecodeOffset;
    result.outputRegister = nodeInfo[outputNode].outputSlotOffset;

    for (MaterialNodeId nodeId : outputList)
    {
        const MaterialNode &node = graph.nodes[nodeId];
        NodeInfo &info = nodeInfo[nodeId];
    }

    return result;
}

} // namespace ybi
