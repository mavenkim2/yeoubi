#include "scene/shader_compilation.h"
#include "core/assert.h"
#include "core/bit_vector.h"
#include "core/thread_context.h"
#include "render/material_types.h"
#include "render/material_opcodes.h"
#include "scene/material.h"
#include <vector>

namespace ybi
{

// struct CompileNode
// {
//     MaterialOpcode type;
//     std::vector<MaterialInput> inputs;
//     int in = 0;
//     bool constantFolded;
//
//     virtual int GetOutputSlot() = 0;
//     CompileNode() = 0;
// };

struct CompileNode
{
    uint32_t outputSlot;
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

void CompileShader(const MaterialGraph &graph, MaterialNodeId outputNode)
{
    // Count number of incoming and outgoing edges
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

    std::vector<uint8_t> compiledBytes;
    std::vector<CompileNode> compiledNodes;

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

        CompileNode compiledNode;
        compiledNode.outputSlot = freeSlot;
        compiledNodes.push_back(compiledNode);
    }

    for (MaterialNodeId nodeId : outputList)
    {
        const MaterialNode &node = graph.nodes[nodeId];
        printf("%s ", ConvertMaterialOpcodeToString(node.type));
        NodeInfo &info = nodeInfo[nodeId];

        printf("output slot: %i\n", info.outputSlotOffset);
    }
}

} // namespace ybi
