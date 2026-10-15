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

#define VMM_KD_TREE_POSITION_STATISTICS_SCALE float(1 << 18)
#define VMM_KD_TREE_WORK_ITEM_SAMPLE_LIMIT 256
#define VMM_KD_TREE_SPLIT_SAMPLE_COUNT (32u * 1024u)

YBI_GPU void CreateKdTreeNodeItems(PathGuidingKdTreeParams &params)
{
    params.numNextNodeItems = 0;
    for (uint32_t activeNodeItemIndex = 0; activeNodeItemIndex < params.numActiveNodeItems;
         ++activeNodeItemIndex)
    {
        const PathGuidingKdTreeNodeItem &nodeItem =
            params.nodeItems[params.activeNodeItems[activeNodeItemIndex]];
        PathGuidingKdTreeNode &node = params.nodes[nodeItem.nodeIndex];
        if (node.childOffset == -1)
        {
            node.childOffset = int(params.numNodes);
            params.numNodes += 2;
            for (uint32_t childIndex = 0; childIndex < 2; ++childIndex)
            {
                PathGuidingKdTreeNode child{};
                child.boundsMin = node.boundsMin;
                child.boundsMax = node.boundsMax;
                if (childIndex == 0)
                    child.boundsMax[node.splitAxis] = node.splitLocation;
                else
                    child.boundsMin[node.splitAxis] = node.splitLocation;
                params.nodes[node.childOffset + childIndex] = child;
            }
        }

        for (uint32_t childIndex = 0; childIndex < 2; ++childIndex)
        {
            PathGuidingKdTreeNodeItem childItem{};
            childItem.nodeIndex = uint32_t(node.childOffset) + childIndex;
            childItem.sampleOffset =
                nodeItem.sampleOffset + (childIndex == 0 ? 0 : nodeItem.numLeft);
            childItem.sampleCount =
                childIndex == 0 ? nodeItem.numLeft : nodeItem.sampleCount - nodeItem.numLeft;
            childItem.treeDepth = nodeItem.treeDepth + 1;
            PathGuidingKdTreeNode &child = params.nodes[childItem.nodeIndex];
            child.numSamples = childItem.sampleCount;
            if (childItem.sampleCount != 0 &&
                (childItem.sampleCount > VMM_KD_TREE_SPLIT_SAMPLE_COUNT || !child.NeedsSplit()))
            {
                params.nextNodeItems[params.numNextNodeItems++] = childItem;
            }
        }
    }
}

YBI_GPU void CreateKdTreeWorkItems(PathGuidingKdTreeParams &params, PathGuidingKdTreeStage stage)
{
    uint32_t totalWorkItemCount = 0;
    params.numActiveNodeItems = 0;
    for (uint32_t nodeItemIndex = 0; nodeItemIndex < params.numNodeItems; nodeItemIndex++)
    {
        PathGuidingKdTreeNodeItem &nodeItem = params.nodeItems[nodeItemIndex];
        PathGuidingKdTreeNode &node = params.nodes[nodeItem.nodeIndex];
        nodeItem.workItemOffset = totalWorkItemCount;
        nodeItem.workItemCount = 0;
        bool matchesStage = node.MatchesStage(stage, params.generation);
        if (stage == PathGuidingKdTreeStage::PositionStatistics)
        {
            node.numSamples = nodeItem.sampleCount;
            matchesStage = matchesStage && nodeItem.sampleCount > VMM_KD_TREE_SPLIT_SAMPLE_COUNT;
        }
        if (!matchesStage || nodeItem.sampleCount == 0)
            continue;

        params.activeNodeItems[params.numActiveNodeItems++] = nodeItemIndex;
        const uint32_t sampleCount = nodeItem.sampleCount;
        const uint32_t workItemCount = CeilDiv(sampleCount, VMM_KD_TREE_WORK_ITEM_SAMPLE_LIMIT);

        nodeItem.workItemCount = workItemCount;
        totalWorkItemCount += workItemCount;

        uint32_t sampleOffset = 0;
        for (uint32_t workItemIndex = nodeItem.workItemOffset;
             workItemIndex < nodeItem.workItemOffset + nodeItem.workItemCount;
             workItemIndex++)
        {
            PathGuidingKdTreeWorkItem &workItem = params.workItems[workItemIndex];
            workItem.nodeIndex = nodeItem.nodeIndex;
            workItem.nodeItemIndex = nodeItemIndex;
            workItem.sampleOffset = nodeItem.sampleOffset + sampleOffset;
            workItem.sampleCount = Min<uint32_t>(VMM_KD_TREE_WORK_ITEM_SAMPLE_LIMIT,
                                                 nodeItem.sampleCount - sampleOffset);
            workItem.leftOffset = 0;

            sampleOffset += VMM_KD_TREE_WORK_ITEM_SAMPLE_LIMIT;
        }
    }
    params.numWorkItems = totalWorkItemCount;
}

YBI_GPU void InitializeNodePositionStatistics(PathGuidingKdTreeParams params,
                                              uint32_t nodeItemIndex)
{
    const PathGuidingKdTreeNodeItem &nodeItem = params.nodeItems[nodeItemIndex];
    PathGuidingKdTreeNode &node = params.nodes[nodeItem.nodeIndex];
    node.boundsMin = Vec3f(INFINITY);
    node.boundsMax = Vec3f(-INFINITY);
    node.meanX = node.meanY = node.meanZ = 0;
    node.varianceX = node.varianceY = node.varianceZ = 0;
    node.numSamples = 0;
}

template <typename Backend>
YBI_GPU void CalculateNodeBounds(PathGuidingKdTreeParams params, uint32_t workItemIndex)
{
    const PathGuidingKdTreeWorkItem &workItem = params.workItems[workItemIndex];
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

template <typename Backend>
YBI_GPU void CalculateNodePositionStatistics(PathGuidingKdTreeParams params,
                                             uint32_t workItemIndex)
{
    const PathGuidingKdTreeWorkItem &workItem = params.workItems[workItemIndex];
    PathGuidingKdTreeNode &node = params.nodes[workItem.nodeIndex];

    const Vec3f boundsHalfExtent = (node.boundsMax - node.boundsMin) * 0.5f;
    const Vec3f boundsCenter = node.boundsMin + boundsHalfExtent;
    YBI_ASSERT(boundsHalfExtent.MinComponent() >= 0.f);
    const Vec3f invHalfExtent(boundsHalfExtent.x > 0.f ? 1.f / boundsHalfExtent.x : 0.f,
                              boundsHalfExtent.y > 0.f ? 1.f / boundsHalfExtent.y : 0.f,
                              boundsHalfExtent.z > 0.f ? 1.f / boundsHalfExtent.z : 0.f);
    int64_t meanX = 0, meanY = 0, meanZ = 0;
    int64_t varianceX = 0, varianceY = 0, varianceZ = 0;
    uint32_t numSamples = 0;

    for (uint32_t sampleIndex = Backend::GetThreadIndex(); sampleIndex < workItem.sampleCount;
         sampleIndex += Backend::GetThreadStride())
    {
        const SampleData &sample =
            params.sampleData[params.sampleIndices[workItem.sampleOffset + sampleIndex]];
        const Vec3f normalized = (sample.position - boundsCenter) * invHalfExtent;
        const Vec3f p(Clamp(normalized.x, -1.f, 1.f),
                      Clamp(normalized.y, -1.f, 1.f),
                      Clamp(normalized.z, -1.f, 1.f));
        meanX += static_cast<int64_t>(VMM_KD_TREE_POSITION_STATISTICS_SCALE * p.x);
        meanY += static_cast<int64_t>(VMM_KD_TREE_POSITION_STATISTICS_SCALE * p.y);
        meanZ += static_cast<int64_t>(VMM_KD_TREE_POSITION_STATISTICS_SCALE * p.z);
        varianceX += static_cast<int64_t>(VMM_KD_TREE_POSITION_STATISTICS_SCALE * p.x * p.x);
        varianceY += static_cast<int64_t>(VMM_KD_TREE_POSITION_STATISTICS_SCALE * p.y * p.y);
        varianceZ += static_cast<int64_t>(VMM_KD_TREE_POSITION_STATISTICS_SCALE * p.z * p.z);
        ++numSamples;
    }
    Backend::ReduceAdd(meanX, node.meanX);
    Backend::ReduceAdd(meanY, node.meanY);
    Backend::ReduceAdd(meanZ, node.meanZ);
    Backend::ReduceAdd(varianceX, node.varianceX);
    Backend::ReduceAdd(varianceY, node.varianceY);
    Backend::ReduceAdd(varianceZ, node.varianceZ);
    Backend::ReduceAdd(numSamples, node.numSamples);
}

YBI_GPU void SelectNodeSplits(PathGuidingKdTreeParams params, uint32_t nodeItemIndex)
{
    const PathGuidingKdTreeNodeItem &nodeItem = params.nodeItems[nodeItemIndex];
    PathGuidingKdTreeNode &node = params.nodes[nodeItem.nodeIndex];
    const Vec3f boundsExtent = node.boundsMax - node.boundsMin;
    if (boundsExtent.MaxComponent() == 0.f)
        return;
    const Vec3f halfExtent = boundsExtent * 0.5f;
    const Vec3f center = node.boundsMin + halfExtent;
    const float invNumSamples = 1.f / float(node.numSamples);
    const Vec3f normalizedMean = Vec3f(float(node.meanX), float(node.meanY), float(node.meanZ)) *
                                 invNumSamples / VMM_KD_TREE_POSITION_STATISTICS_SCALE;
    Vec3f mean = normalizedMean * halfExtent + center;
    const Vec3f halfBinSize = halfExtent * (0.5f / VMM_KD_TREE_POSITION_STATISTICS_SCALE);
    for (uint32_t axis = 0; axis < 3; ++axis)
    {
        const float inset = halfBinSize[axis];
        if (mean[axis] <= node.boundsMin[axis])
            mean[axis] = node.boundsMin[axis] + inset;
        if (mean[axis] >= node.boundsMax[axis])
            mean[axis] = node.boundsMax[axis] - inset;
    }

    Vec3f variance = Vec3f(float(node.varianceX), float(node.varianceY), float(node.varianceZ)) /
                     VMM_KD_TREE_POSITION_STATISTICS_SCALE * invNumSamples;
    variance -= normalizedMean * normalizedMean;
    variance = Max(variance, Vec3f(0.f)) * (halfExtent * halfExtent);
    variance = Min(variance, boundsExtent * boundsExtent);

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
    if (mean[splitAxis] <= node.boundsMin[splitAxis])
        mean[splitAxis] = nextafterf(node.boundsMin[splitAxis], node.boundsMax[splitAxis]);
    node.splitAxis = int(splitAxis);
    node.splitLocation = mean[splitAxis];
}

template <typename Backend>
YBI_GPU void PartitionCount(PathGuidingKdTreeParams params, uint32_t workItemIndex)
{
    PathGuidingKdTreeWorkItem &workItem = params.workItems[workItemIndex];
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
    Backend::ReduceAdd(numLefts, workItem.leftOffset);
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
            count = params.workItems[nodeItem.workItemOffset + workItemOffset].leftOffset;
        }

        const auto scan = Backend::ExclusiveScanAdd(count);
        if (valid)
        {
            params.workItems[nodeItem.workItemOffset + workItemOffset].leftOffset =
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
    const PathGuidingKdTreeWorkItem &workItem = params.workItems[workItemIndex];
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
            const uint32_t outputOffset = goesLeft ? workItem.leftOffset + localLeft
                                                   : nodeItem.numLeft +
                                                         (precedingSamples - workItem.leftOffset) +
                                                         (sampleOffset - localLeft);
            newSampleIndices[nodeItem.sampleOffset + outputOffset] = sampleIndex;
        }
        previousLeft += scan.total;
    }
    if (threadIndex == 0 && workItemIndex == nodeItem.workItemOffset)
        node.partitionGeneration = params.generation;
}

} // namespace ybi
