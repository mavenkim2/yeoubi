#pragma once

#include "kernel/guiding/types.h"

namespace ybi
{

YBI_CPU_GPU const VMMDistribution *TraverseKdTree(const PathGuidingKdTreeParams &params,
                                                  const Vec3f &position)
{
    if (!params.nodes || params.numNodes == 0 || !params.vmmDistributions)
        return nullptr;

    const PathGuidingKdTreeNode &root = params.nodes[0];
    for (uint32_t axis = 0; axis < 3; ++axis)
    {
        if (!(position[axis] >= root.cellBoundsMin[axis] &&
              position[axis] <= root.cellBoundsMax[axis]))
            return nullptr;
    }

    uint32_t nodeIndex = 0;
    while (!params.nodes[nodeIndex].IsGuidingCell())
    {
        const PathGuidingKdTreeNode &node = params.nodes[nodeIndex];
        if (!node.HasChildren())
            return nullptr;

        YBI_ASSERT(node.splitAxis >= 0 && node.splitAxis < 3);
        YBI_ASSERT(uint32_t(node.childOffset) + 1 < params.numNodes);
        nodeIndex =
            uint32_t(node.childOffset) + (position[node.splitAxis] < node.splitLocation ? 0u : 1u);
    }

    const uint32_t distributionIndex = uint32_t(params.nodes[nodeIndex].guidingCellIndex);
    YBI_ASSERT(distributionIndex < params.numVMMs);
    const VMMDistribution &distribution = params.vmmDistributions[distributionIndex];
    return distribution.initialized && distribution.numComponents > 0 ? &distribution : nullptr;
}

} // namespace ybi
