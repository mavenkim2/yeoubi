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
                child.boundsMin = node.boundsMin;
                child.boundsMax = node.boundsMax;
                if (childIndex == 0)
                    child.boundsMax[node.splitAxis] = node.splitLocation;
                else
                    child.boundsMin[node.splitAxis] = node.splitLocation;
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
            PathGuidingKdTreeNode &child = params.nodes[childItem.nodeIndex];
            child.numSamples = childItem.sampleCount;
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
        PathGuidingKdTreeNode &node = params.nodes[nodeItem.nodeIndex];
        nodeItem.workItemOffset = totalWorkItemCount;
        nodeItem.workItemCount = 0;
        const bool matchesStage = node.MatchesStage(stage, params.generation);
        if (stage == PathGuidingKdTreeStage::PositionStatistics)
            node.numSamples = nodeItem.sampleCount;
        if (!matchesStage || nodeItem.sampleCount == 0)
            continue;

        if (stage == PathGuidingKdTreeStage::PositionStatistics)
        {
            node.boundsMin = Vec3f(INFINITY);
            node.boundsMax = Vec3f(-INFINITY);
        }
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
    Vec3f boundsMin(INFINITY);
    Vec3f boundsMax(-INFINITY);
    for (uint32_t sampleOffset = Backend::GetThreadIndex(); sampleOffset < workItem.sampleCount;
         sampleOffset += Backend::GetThreadStride())
    {
        const uint32_t sampleIndex = params.sampleIndices[workItem.sampleOffset + sampleOffset];
        const Vec3f &position = params.sampleData[sampleIndex].position;
        boundsMin = Min(boundsMin, position);
        boundsMax = Max(boundsMax, position);
    }
    Backend::ReduceMin(boundsMin, node.boundsMin);
    Backend::ReduceMax(boundsMax, node.boundsMax);
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
        const SampleData &sample =
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

YBI_GPU void StoreNodePositionStatistics(PathGuidingKdTreeParams params,
                                         uint32_t nodeItemIndex,
                                         uint32_t workItemIndex,
                                         const PathGuidingKdTreePositionStatistics &statistics,
                                         PathGuidingKdTreePositionStatistics *output)
{
    output[workItemIndex] = statistics;
    const PathGuidingKdTreeNodeItem &nodeItem = params.nodeItems[nodeItemIndex];
    if (nodeItem.workItemCount == 1)
    {
        PathGuidingKdTreeNode &node = params.nodes[nodeItem.nodeIndex];
        node.meanX = statistics.mean.x;
        node.meanY = statistics.mean.y;
        node.meanZ = statistics.mean.z;
        node.varianceX = statistics.variance.x;
        node.varianceY = statistics.variance.y;
        node.varianceZ = statistics.variance.z;
        node.numSamples = statistics.numSamples;
    }
}

YBI_GPU void SelectNodeSplits(PathGuidingKdTreeParams params, uint32_t nodeItemIndex)
{
    const PathGuidingKdTreeNodeItem &nodeItem = params.nodeItems[nodeItemIndex];
    PathGuidingKdTreeNode &node = params.nodes[nodeItem.nodeIndex];
    if (params.splitDecision == PathGuidingSplitDecision::SampleCount)
    {
        if (node.numSamples <= VMM_KD_TREE_SPLIT_SAMPLE_COUNT)
            return;
    }
    else if (node.numSamples < VMM_KD_TREE_LOOKAHEAD_SPLIT_SAMPLE_COUNT ||
             GetLookaheadDepth(params, node) >= params.maxLookaheadDepth)
    {
        return;
    }
    const Vec3f boundsExtent = node.boundsMax - node.boundsMin;
    if (boundsExtent.MaxComponent() == 0.f)
        return;
    Vec3f mean(node.meanX, node.meanY, node.meanZ);
    const Vec3f variance(node.varianceX, node.varianceY, node.varianceZ);

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
    mean[splitAxis] = Min(mean[splitAxis], node.boundsMax[splitAxis]);
    if (mean[splitAxis] <= node.boundsMin[splitAxis])
        mean[splitAxis] = nextafterf(node.boundsMin[splitAxis], node.boundsMax[splitAxis]);
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
        const SampleData &sample =
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
