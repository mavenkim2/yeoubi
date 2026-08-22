#include "scene/shader_compilation.h"
#include "core/assert.h"
#include "core/bit_vector.h"
#include "core/thread_context.h"
#include "render/material_types.h"
#include "render/material_opcodes.h"
#include "scene/material.h"
#include <unordered_map>
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
static void IterateMaterialNodeInputs(MaterialNode *node, Func &&func)
{
    for (const MaterialInput &input : node->inputs)
    {
        std::visit(
            [&](const auto &value) {
                using T = std::decay_t<decltype(value)>;
                if constexpr (std::is_same_v<T, std::shared_ptr<MaterialNode>>)
                {
                    func(value.get());
                }
            },
            input.value);
    }
}

void CompileShader(std::shared_ptr<MaterialNode> outputNode)
{
    // Count number of incoming and outgoing edges
    struct NodeInfo
    {
        std::vector<MaterialNode *> outNodes;
        int in = 0;
        int out = 0;

        // Used when compiling
        int inputSlotOffset;
        int numInputSlots;

        int outputSlotOffset;
        int numOutputSlots;
    };

    YBI_ASSERT(outputNode);

    std::unordered_map<MaterialNode *, NodeInfo> nodeMap;
    std::vector<MaterialNode *> discoveryOrder;
    std::vector<MaterialNode *> stack;

    stack.push_back(outputNode.get());
    nodeMap.try_emplace(outputNode.get());
    discoveryOrder.push_back(outputNode.get());

    while (!stack.empty())
    {
        MaterialNode *node = stack.back();
        stack.pop_back();
        IterateMaterialNodeInputs(node, [&](MaterialNode *parent) {
            nodeMap[node].in++;
            auto [parentEntry, inserted] = nodeMap.try_emplace(parent);
            parentEntry->second.outNodes.push_back(node);
            parentEntry->second.out++;

            if (inserted)
            {
                discoveryOrder.push_back(parent);
                stack.push_back(parent);
            }
        });
    }

    // Topological sort the nodes
    std::vector<MaterialNode *> queue;
    for (MaterialNode *node : discoveryOrder)
    {
        if (nodeMap.at(node).in == 0)
        {
            queue.push_back(node);
        }
    }

    std::vector<MaterialNode *> outputList;
    size_t index = 0;
    while (index != queue.size())
    {
        MaterialNode *node = queue[index++];
        outputList.push_back(node);
        auto entry = nodeMap.find(node);
        YBI_ASSERT(entry != nodeMap.end());
        for (MaterialNode *neighbor : entry->second.outNodes)
        {
            auto neighborEntry = nodeMap.find(neighbor);
            YBI_ASSERT(neighborEntry != nodeMap.end());
            YBI_ASSERT(neighborEntry->second.in);
            if (--neighborEntry->second.in == 0)
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

    for (MaterialNode *node : outputList)
    {
        size_t numWords = GetMaterialOpcodeOutputRegisterCount(node->type);

#if 0
        // See if node computes a constant
        bool constant = true;
        for (MaterialInput &input : node->inputs)
        {
            std::visit(
                [&](const auto &value) {
                    using T = std::decay_t<decltype(value)>;

                    if constexpr (std::is_same_v<T, std::shared_ptr<MaterialNode>>)
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
            printf("%s\n", ConvertMaterialOpcodeToString(node->type));
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
                  ConvertMaterialOpcodeToString(outputNode->type));

        // Free slots from any parents whose output is no longer used
        IterateMaterialNodeInputs(node, [&](MaterialNode *parent) {
            NodeInfo &info = nodeMap[parent];
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

        NodeInfo &info = nodeMap[node];
        info.outputSlotOffset = freeSlot;
        info.numOutputSlots = numWords;

        CompileNode node;
        node.outputSlot = freeSlot;
        compiledNodes.push_back(node);
    }

    for (MaterialNode *node : outputList)
    {
        printf("%s ", ConvertMaterialOpcodeToString(node->type));
        NodeInfo &info = nodeMap[node];

        printf("output slot: %i\n", info.outputSlotOffset);
    }
}

} // namespace ybi
