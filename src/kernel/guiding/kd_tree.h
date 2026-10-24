// Implementation adapted from https://github.com/OpenPathGuidingLibrary/openpgl.
// Copyright 2020 Intel Corporation
// SPDX-License-Identifier: Apache-2.0
// Full license: licenses/OpenPGL.txt (relative to the repository root).
#pragma once

#include "core/math.h"
#include "core/vec3f.h"
#include "device/cuda_compat.h"
#include "kernel/guiding/types.h"

namespace ybi
{

#define VMM_KD_TREE_WORK_ITEM_SAMPLE_LIMIT 256
#define VMM_KD_TREE_SPLIT_SAMPLE_COUNT (32u * 1024u)
#define VMM_KD_TREE_LOOKAHEAD_SPLIT_SAMPLE_COUNT 1000u

YBI_GPU void CollectVMMUpdateNodes(PathGuidingKdTreeParams &params)
{
    params.numVMMUpdateNodes = 0;
    for (uint32_t nodeIndex = 0; nodeIndex < params.numNodes; ++nodeIndex)
    {
        if (params.nodes[nodeIndex].IsGuidingCell() && params.nodeItems[nodeIndex].sampleCount > 0)
            params.vmmUpdateNodeIndices[params.numVMMUpdateNodes++] = nodeIndex;
    }
}

YBI_GPU void CopySampleCountDistributions(PathGuidingKdTreeParams &params)
{
    for (uint32_t frontierIndex = 0; frontierIndex < params.numNodeItemIndices; ++frontierIndex)
    {
        PathGuidingKdTreeNode &node = params.nodes[params.nodeItemIndices[frontierIndex]];
        if (!node.IsGuidingCell() || !node.HasChildren())
            continue;

        PathGuidingKdTreeNode &left = params.nodes[node.childOffset];
        PathGuidingKdTreeNode &right = params.nodes[node.childOffset + 1];
        left.guidingCellIndex = node.guidingCellIndex;
        right.guidingCellIndex = int(params.numVMMs++);
        params.vmmDistributions[right.guidingCellIndex] =
            params.vmmDistributions[left.guidingCellIndex];
        params.vmmDistributions[left.guidingCellIndex].split = true;
        params.vmmDistributions[right.guidingCellIndex].split = true;
        if (params.vmmStatistics)
            params.vmmStatistics[right.guidingCellIndex] =
                params.vmmStatistics[left.guidingCellIndex];
        if (params.vmmSplitStatistics)
            params.vmmSplitStatistics[right.guidingCellIndex] =
                params.vmmSplitStatistics[left.guidingCellIndex];
        node.guidingCellIndex = -1;
    }
}

YBI_GPU uint32_t GetLookaheadDepth(const PathGuidingKdTreeParams &params,
                                   const PathGuidingKdTreeNode &node)
{
    if (node.IsGuidingCell())
        return 0;
    YBI_ASSERT(node.IsLookahead());
    return node.treeDepth - params.nodes[node.ancestorGuidingCellIndex].treeDepth;
}

YBI_GPU bool IsGuidingCellMarkedForSplit(const PathGuidingKdTreeParams &params, uint32_t nodeIndex)
{
    return (params.guidingCellSplitMask[nodeIndex >> 5] & (1u << (nodeIndex & 31u))) != 0;
}

YBI_GPU bool LookaheadContinuationCriteria(const PathGuidingKdTreeParams &params,
                                           uint32_t nodeItemIndex)
{
    const PathGuidingKdTreeNodeItem &nodeItem = params.nodeItems[nodeItemIndex];
    const PathGuidingKdTreeNode &node = params.nodes[nodeItem.nodeIndex];
    if (!node.IsGuidingCell() && !node.IsLookahead())
        return true;

    const uint32_t ancestorIndex =
        node.IsGuidingCell() ? nodeItem.nodeIndex : uint32_t(node.ancestorGuidingCellIndex);
    // If the ancestor needs to be split, then its children will be considered in the
    // next pass instead of the current node's children.
    return params.nodes[ancestorIndex].IsGuidingCell() &&
           !IsGuidingCellMarkedForSplit(params, ancestorIndex);
}

YBI_GPU void CreateKdTreeNodes(PathGuidingKdTreeParams &params,
                               PathGuidingSplitDecision splitDecision)
{
    for (uint32_t frontierIndex = 0; frontierIndex < params.numNodeItemIndices; ++frontierIndex)
    {
        const PathGuidingKdTreeNodeItem &nodeItem =
            params.nodeItems[params.nodeItemIndices[frontierIndex]];
        PathGuidingKdTreeNode &node = params.nodes[nodeItem.nodeIndex];
        if (nodeItem.sampleCount == 0 || node.NeedsSplit())
            continue;
        if (!node.HasChildren())
        {
            node.childOffset = int(params.numNodes);
            params.numNodes += 2;
            for (uint32_t childIndex = 0; childIndex < 2; ++childIndex)
            {
                PathGuidingKdTreeNode child{};
                child.treeDepth = node.treeDepth + 1;
                if (splitDecision == PathGuidingSplitDecision::Lookahead)
                {
                    child.ancestorGuidingCellIndex = node.IsGuidingCell()
                                                         ? int(nodeItem.nodeIndex)
                                                         : node.ancestorGuidingCellIndex;
                }
                child.cellBoundsMin = node.cellBoundsMin;
                child.cellBoundsMax = node.cellBoundsMax;
                if (childIndex == 0)
                    child.cellBoundsMax[node.splitAxis] = node.splitLocation;
                else
                    child.cellBoundsMin[node.splitAxis] = node.splitLocation;
                params.nodes[node.childOffset + childIndex] = child;
                params.nodeItems[node.childOffset + childIndex] = {};
            }
        }

        for (uint32_t childIndex = 0; childIndex < 2; ++childIndex)
        {
            const uint32_t childNodeIndex = uint32_t(node.childOffset) + childIndex;
            PathGuidingKdTreeNodeItem &childItem = params.nodeItems[childNodeIndex];
            childItem.nodeIndex = childNodeIndex;
            childItem.sampleOffset =
                nodeItem.sampleOffset + (childIndex == 0 ? 0 : nodeItem.numLeft);
            childItem.sampleCount =
                childIndex == 0 ? nodeItem.numLeft : nodeItem.sampleCount - nodeItem.numLeft;
        }
    }
}

YBI_GPU void CreateKdTreeNodeItems(PathGuidingKdTreeParams &params,
                                   PathGuidingSplitDecision splitDecision)
{
    params.numNextNodeItemIndices = 0;
    for (uint32_t frontierIndex = 0; frontierIndex < params.numNodeItemIndices; ++frontierIndex)
    {
        const uint32_t nodeItemIndex = params.nodeItemIndices[frontierIndex];
        const PathGuidingKdTreeNodeItem &nodeItem = params.nodeItems[nodeItemIndex];
        const PathGuidingKdTreeNode &node = params.nodes[nodeItem.nodeIndex];
        if (nodeItem.sampleCount == 0 || !node.HasChildren())
            continue;
        if (splitDecision == PathGuidingSplitDecision::Lookahead &&
            !LookaheadContinuationCriteria(params, nodeItemIndex))
            continue;
        for (uint32_t childIndex = 0; childIndex < 2; ++childIndex)
        {
            const uint32_t childNodeIndex = uint32_t(node.childOffset) + childIndex;
            if (params.nodeItems[childNodeIndex].sampleCount != 0)
            {
                params.nextNodeItemIndices[params.numNextNodeItemIndices++] = childNodeIndex;
            }
        }
    }
}

YBI_GPU void CreateKdTreeWorkItems(PathGuidingKdTreeParams &params, PathGuidingKdTreeStage stage)
{
    uint32_t totalWorkItemCount = 0;
    params.numActiveNodeItems = 0;
    params.maxWorkItemCount = 0;
    for (uint32_t frontierIndex = 0; frontierIndex < params.numNodeItemIndices; ++frontierIndex)
    {
        const uint32_t nodeItemIndex = params.nodeItemIndices[frontierIndex];
        PathGuidingKdTreeNodeItem &nodeItem = params.nodeItems[nodeItemIndex];
        const PathGuidingKdTreeNode &node = params.nodes[nodeItem.nodeIndex];
        nodeItem.workItemOffset = totalWorkItemCount;
        nodeItem.workItemCount = 0;
        const bool matchesStage = node.MatchesStage(stage, params.generation);
        if (!matchesStage || nodeItem.sampleCount == 0)
            continue;
        params.activeNodeItems[params.numActiveNodeItems++] = nodeItemIndex;
        const uint32_t sampleCount = nodeItem.sampleCount;
        const uint32_t workItemCount = CeilDiv(sampleCount, VMM_KD_TREE_WORK_ITEM_SAMPLE_LIMIT);

        nodeItem.workItemCount = workItemCount;
        params.maxWorkItemCount = Max(params.maxWorkItemCount, workItemCount);
        totalWorkItemCount += workItemCount;
    }
    params.numWorkItems = totalWorkItemCount;
}

YBI_GPU uint32_t FindKdTreeNodeItemIndex(PathGuidingKdTreeParams params, uint32_t workItemIndex)
{
    uint32_t begin = 0;
    uint32_t end = params.numActiveNodeItems;
    while (begin + 1 < end)
    {
        const uint32_t middle = begin + (end - begin) / 2;
        const PathGuidingKdTreeNodeItem &nodeItem =
            params.nodeItems[params.activeNodeItems[middle]];
        if (workItemIndex < nodeItem.workItemOffset)
            end = middle;
        else
            begin = middle;
    }

    return params.activeNodeItems[begin];
}

YBI_GPU PathGuidingKdTreeWorkItem GetKdTreeWorkItem(PathGuidingKdTreeParams params,
                                                    uint32_t workItemIndex)
{
    const uint32_t nodeItemIndex = FindKdTreeNodeItemIndex(params, workItemIndex);
    const PathGuidingKdTreeNodeItem &nodeItem = params.nodeItems[nodeItemIndex];
    const uint32_t sampleOffset =
        (workItemIndex - nodeItem.workItemOffset) * VMM_KD_TREE_WORK_ITEM_SAMPLE_LIMIT;
    PathGuidingKdTreeWorkItem workItem{};
    workItem.nodeIndex = nodeItem.nodeIndex;
    workItem.nodeItemIndex = nodeItemIndex;
    workItem.sampleOffset = nodeItem.sampleOffset + sampleOffset;
    workItem.sampleCount =
        Min<uint32_t>(VMM_KD_TREE_WORK_ITEM_SAMPLE_LIMIT, nodeItem.sampleCount - sampleOffset);
    return workItem;
}

template <typename Backend>
YBI_GPU void CalculateNodeBounds(PathGuidingKdTreeParams params, uint32_t workItemIndex)
{
    const PathGuidingKdTreeWorkItem workItem = GetKdTreeWorkItem(params, workItemIndex);
    PathGuidingKdTreeNode &node = params.nodes[workItem.nodeIndex];
    Vec3f boundsMin(FLOAT_INFINITY);
    Vec3f boundsMax(-FLOAT_INFINITY);
    for (uint32_t sampleOffset = Backend::GetThreadIndex(); sampleOffset < workItem.sampleCount;
         sampleOffset += Backend::GetThreadStride())
    {
        const uint32_t sampleIndex = params.sampleIndices[workItem.sampleOffset + sampleOffset];
        const Vec3f &position = params.sampleData[sampleIndex].position;
        boundsMin = Min(boundsMin, position);
        boundsMax = Max(boundsMax, position);
    }
    Backend::ReduceMin(boundsMin, node.positionBoundsMin);
    Backend::ReduceMax(boundsMax, node.positionBoundsMax);
}

YBI_GPU PathGuidingKdTreePositionStatistics
CalculateNodePositionStatistics(PathGuidingKdTreeParams params,
                                const PathGuidingKdTreeWorkItem &workItem,
                                uint32_t threadIndex)
{
    PathGuidingKdTreePositionStatistics result = {};
    for (uint32_t sampleOffset = threadIndex; sampleOffset < workItem.sampleCount;
         sampleOffset += VMM_BLOCK_SIZE)
    {
        const PathGuidingTrainingSample &sample =
            params.sampleData[params.sampleIndices[workItem.sampleOffset + sampleOffset]];
        const PathGuidingKdTreePositionStatistics value{sample.position, Vec3f(0.f), 1};
        result = result + value;
    }
    return result;
}

YBI_GPU void CreateKdTreeReductionWorkItems(PathGuidingKdTreeParams &params)
{
    uint32_t totalWorkItemCount = 0;
    for (uint32_t activeIndex = 0; activeIndex < params.numActiveNodeItems; ++activeIndex)
    {
        PathGuidingKdTreeNodeItem &nodeItem =
            params.nodeItems[params.activeNodeItems[activeIndex]];
        nodeItem.reductionInputOffset = nodeItem.workItemOffset;
        nodeItem.reductionInputCount = nodeItem.workItemCount;
        nodeItem.workItemOffset = totalWorkItemCount;
        nodeItem.workItemCount = CeilDiv(nodeItem.reductionInputCount, VMM_BLOCK_SIZE);
        totalWorkItemCount += nodeItem.workItemCount;
    }
    params.numWorkItems = totalWorkItemCount;
}

template <typename T>
YBI_GPU T LoadKdTreeReductionValue(const PathGuidingKdTreeNodeItem &nodeItem,
                                   uint32_t workItemIndex,
                                   uint32_t threadIndex,
                                   const T *input)
{
    const uint32_t offset =
        (workItemIndex - nodeItem.workItemOffset) * VMM_BLOCK_SIZE + threadIndex;
    if (offset >= nodeItem.reductionInputCount)
        return {};
    return input[nodeItem.reductionInputOffset + offset];
}

YBI_GPU void AccumulateNodePositionStatistics(const PathGuidingKdTreeParams &params,
                                              uint32_t nodeItemIndex,
                                              const PathGuidingKdTreePositionStatistics *input)
{
    PathGuidingKdTreeNodeItem &nodeItem = params.nodeItems[nodeItemIndex];
    // TODO: Avoid redundant bounds/statistics reductions on revisits while still
    // allowing split selection to reconsider the node after promotion.
    if (nodeItem.positionStatisticsAccumulated)
        return;

    // Only commit after every reduction level has finished. Nodes with one work
    // item may otherwise appear repeatedly while larger nodes are still reduced.
    YBI_ASSERT(nodeItem.workItemCount == 1);
    const auto &statistics = input[nodeItem.workItemOffset];
    YBI_ASSERT(statistics.numSamples == nodeItem.sampleCount);
    PathGuidingKdTreeNode &node = params.nodes[nodeItem.nodeIndex];
    node.positionStatistics = node.positionStatistics + statistics;
    nodeItem.positionStatisticsAccumulated = true;
}

YBI_GPU void SelectNodeSplits(PathGuidingKdTreeParams params, uint32_t nodeItemIndex)
{
    const PathGuidingKdTreeNodeItem &nodeItem = params.nodeItems[nodeItemIndex];
    PathGuidingKdTreeNode &node = params.nodes[nodeItem.nodeIndex];
    const auto &statistics = node.positionStatistics;
    if (params.splitDecision == PathGuidingSplitDecision::SampleCount)
    {
        if (statistics.numSamples <= VMM_KD_TREE_SPLIT_SAMPLE_COUNT)
            return;
    }
    else if (statistics.numSamples < VMM_KD_TREE_LOOKAHEAD_SPLIT_SAMPLE_COUNT ||
             GetLookaheadDepth(params, node) >= params.maxLookaheadDepth)
    {
        return;
    }
    const Vec3f boundsExtent = node.positionBoundsMax - node.positionBoundsMin;
    if (boundsExtent.MaxComponent() == 0.f)
        return;
    Vec3f mean = statistics.mean;
    const Vec3f variance = statistics.variance;

    const bool equalVariance = variance.x == variance.y && variance.y == variance.z;
    uint32_t splitAxis;
    if (equalVariance && boundsExtent.x == boundsExtent.y && boundsExtent.y == boundsExtent.z)
    {
        splitAxis = 0;
    }
    else
    {
        const Vec3f axisValues = equalVariance ? boundsExtent : variance;
        const uint32_t xyAxis = axisValues.y > axisValues.x ? 1 : 0;
        splitAxis = axisValues[xyAxis] > axisValues.z ? xyAxis : 2;
    }
    mean[splitAxis] = Min(mean[splitAxis], node.positionBoundsMax[splitAxis]);
    if (mean[splitAxis] <= node.positionBoundsMin[splitAxis])
        mean[splitAxis] =
            nextafterf(node.positionBoundsMin[splitAxis], node.positionBoundsMax[splitAxis]);
    node.splitAxis = int(splitAxis);
    node.splitLocation = mean[splitAxis];
}

template <typename Backend>
YBI_GPU void PartitionCount(PathGuidingKdTreeParams params, uint32_t workItemIndex)
{
    const PathGuidingKdTreeWorkItem workItem = GetKdTreeWorkItem(params, workItemIndex);
    const PathGuidingKdTreeNode &node = params.nodes[workItem.nodeIndex];

    uint32_t numLefts = 0;
    for (uint32_t sampleIndex = Backend::GetThreadIndex(); sampleIndex < workItem.sampleCount;
         sampleIndex += Backend::GetThreadStride())
    {
        const PathGuidingTrainingSample &sample =
            params.sampleData[params.sampleIndices[workItem.sampleOffset + sampleIndex]];
        const bool predicate = sample.position[node.splitAxis] < node.splitLocation;
        numLefts += predicate;
    }
    const auto scan = Backend::ExclusiveScanAdd(numLefts);
    if (Backend::GetThreadIndex() == 0)
        params.workItemLeftOffsets[workItemIndex] = scan.total;
}

template <typename Backend>
YBI_GPU void PartitionExclusiveScan(PathGuidingKdTreeParams params, uint32_t nodeItemIndex)
{
    PathGuidingKdTreeNodeItem &nodeItem = params.nodeItems[nodeItemIndex];
    const uint32_t threadIndex = Backend::GetThreadIndex();
    const uint32_t threadStride = Backend::GetThreadStride();

    uint32_t offset = 0;
    for (uint32_t tileOffset = 0; tileOffset < nodeItem.workItemCount; tileOffset += threadStride)
    {
        const uint32_t workItemOffset = tileOffset + threadIndex;
        const bool valid = workItemOffset < nodeItem.workItemCount;
        uint32_t count = 0;
        if (valid)
        {
            count = params.workItemLeftOffsets[nodeItem.workItemOffset + workItemOffset];
        }

        const auto scan = Backend::ExclusiveScanAdd(count);
        if (valid)
        {
            params.workItemLeftOffsets[nodeItem.workItemOffset + workItemOffset] =
                offset + scan.prefix;
        }
        offset += scan.total;
    }
    if (threadIndex == 0)
        nodeItem.numLeft = offset;
}

template <typename Backend>
YBI_GPU void PartitionSamples(PathGuidingKdTreeParams params,
                              uint32_t workItemIndex,
                              uint32_t *newSampleIndices)
{
    const uint32_t threadIndex = Backend::GetThreadIndex();
    const uint32_t threadStride = Backend::GetThreadStride();
    const PathGuidingKdTreeWorkItem workItem = GetKdTreeWorkItem(params, workItemIndex);
    const uint32_t leftOffset = params.workItemLeftOffsets[workItemIndex];
    PathGuidingKdTreeNode &node = params.nodes[workItem.nodeIndex];
    const PathGuidingKdTreeNodeItem &nodeItem = params.nodeItems[workItem.nodeItemIndex];
    const uint32_t precedingSamples = workItem.sampleOffset - nodeItem.sampleOffset;
    uint32_t previousLeft = 0;

    for (uint32_t tileOffset = 0; tileOffset < workItem.sampleCount; tileOffset += threadStride)
    {
        const uint32_t sampleOffset = tileOffset + threadIndex;
        const bool valid = sampleOffset < workItem.sampleCount;
        uint32_t sampleIndex = 0;
        bool goesLeft = false;
        if (valid)
        {
            sampleIndex = params.sampleIndices[workItem.sampleOffset + sampleOffset];
            goesLeft =
                params.sampleData[sampleIndex].position[node.splitAxis] < node.splitLocation;
        }

        const auto scan = Backend::ExclusiveScanAdd(uint32_t(goesLeft));
        const uint32_t localLeft = previousLeft + scan.prefix;
        if (valid)
        {
            const uint32_t outputOffset = goesLeft ? leftOffset + localLeft
                                                   : nodeItem.numLeft +
                                                         (precedingSamples - leftOffset) +
                                                         (sampleOffset - localLeft);
            newSampleIndices[nodeItem.sampleOffset + outputOffset] = sampleIndex;
        }
        previousLeft += scan.total;
    }
    if (threadIndex == 0 && workItemIndex == nodeItem.workItemOffset)
        node.partitionGeneration = params.generation;
}

} // namespace ybi
