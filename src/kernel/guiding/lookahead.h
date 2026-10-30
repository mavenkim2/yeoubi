// Implementation adapted from https://github.com/JamesZFS/OpenPGL-spatial-subdiv.
// Copyright 2020 Intel Corporation
// SPDX-License-Identifier: Apache-2.0
// Full license: licenses/OpenPGL.txt (relative to the repository root).
#pragma once

#include "core/atomic.h"
#include "core/octahedral.h"
#include "core/vec3f.h"
#include "device/cuda_compat.h"
#include "kernel/guiding/kd_tree.h"
#include "kernel/guiding/types.h"
#include "kernel/guiding/vmf.h"
#include "kernel/texture/sampling.h"

namespace ybi
{

// Fengshi Zheng, Christoph Peters, Sebastian Herholz, Marco Manzi, and Elmar Eisemann.
// "Illumination-Aware Spatial Subdivision for Path Guiding", EGSR, 2026.
// https://doi.org/10.2312/sr.20261000

#define VMM_KD_TREE_CRITICAL_Z_SCORE 3.7190164855f
#define VMM_KD_TREE_MIN_NUM_SAMPLES 1000u

YBI_GPU Vec3f CalculateMeanDirection(const LookaheadCell &cell)
{
    return Normalize(cell.radianceWeightedDirectionSum);
}

YBI_GPU float CalculateMeanDirectionInvKappa(const LookaheadCell &cell)
{
    const Vec3f r = cell.radianceWeightedDirectionSum / cell.radianceFirstMoment;
    const float meanCosine = Min(Length(r), 1.f);
    if (meanCosine == 1.f)
        return 0.f; // Infinite concentration.

    const float kappa = MeanCosineToKappa(meanCosine);
    const float kappaCothMinusOne = kappa * meanCosine;
    const float yMinusOne = (cell.radianceFirstMoment / cell.radianceSecondMoment) *
                            cell.radianceFirstMoment * kappaCothMinusOne;
    const float invY = 1.f / (1.f + yMinusOne);
    // Eq. (4), evaluated in 1/y to avoid polynomial overflow.
    const float numerator = 1.f + invY * (1.69934861f + invY * (5.38753272f + invY * 9.85021305f));
    const float denominator = 1.f + invY * (0.67453491f + invY * 4.31180006f);
    return Sqrt((invY / yMinusOne) * (denominator / numerator));
}

YBI_GPU float CalculateSampleMean(const LookaheadCell &cell)
{
    return cell.radianceFirstMoment / cell.numSamples;
}

YBI_GPU float CalculateSampleVariance(const LookaheadCell &cell)
{
    // Variance of the estimated mean, not of an individual sample.
    const double n = cell.numSamples;
    const double b1 = cell.radianceFirstMoment;
    const double numerator = n * cell.radianceSecondMoment - b1 * b1;
    // Accumulated float moments can produce a small negative value through rounding.
    return static_cast<float>(Max(0.0, numerator) / (n * n * n));
}

YBI_GPU LookaheadCell CalculateDelta(const LookaheadCell &lhs, const LookaheadCell &rhs)
{
    YBI_ASSERT(rhs.numSamples <= lhs.numSamples);
    LookaheadCell result = {lhs.radianceFirstMoment - rhs.radianceFirstMoment,
                            lhs.radianceSecondMoment - rhs.radianceSecondMoment,
                            lhs.numSamples - rhs.numSamples,
                            lhs.radianceWeightedDirectionSum - rhs.radianceWeightedDirectionSum};
#if defined(YBI_ENABLE_GUIDING_RGB_STATS)
    result.rgbRadianceFirstMoment = lhs.rgbRadianceFirstMoment - rhs.rgbRadianceFirstMoment;
#endif

    return result;
}

YBI_GPU PathGuidingKdTreeLookaheadStatistics
UpdateLookaheadCell(const PathGuidingKdTreeParams &params,
                    const PathGuidingKdTreeWorkItem &workItem,
                    uint32_t threadIndex)
{
    PathGuidingKdTreeLookaheadStatistics result = {};
    for (uint32_t sampleOffset = threadIndex; sampleOffset < workItem.sampleCount;
         sampleOffset += VMM_BLOCK_SIZE)
    {
        const uint32_t sampleIndex = params.sampleIndices[workItem.sampleOffset + sampleOffset];
        const PathGuidingTrainingSample &sample = params.sampleData[sampleIndex];
        YBI_ASSERT(sample.pdf > 0.f);
        const float Li = sample.GetIncidentRadiance() / sample.pdf;
        const Vec3f direction = DecodeOctahedral(sample.octahedralDirection);
        PathGuidingKdTreeLookaheadStatistics value{Li, Sqr(Li), direction * Li, 1};
#if defined(YBI_ENABLE_GUIDING_RGB_STATS)
        value.meanRadianceRGB = sample.incidentRadiance / sample.pdf;
#endif
        result = result + value;
    }
    return result;
}

YBI_GPU void AccumulateLookaheadCell(const PathGuidingKdTreeParams &params,
                                     uint32_t nodeItemIndex,
                                     const PathGuidingKdTreeLookaheadStatistics *input)
{
    const PathGuidingKdTreeNodeItem &nodeItem = params.nodeItems[nodeItemIndex];
    const PathGuidingKdTreeLookaheadStatistics &statistics = input[nodeItem.workItemOffset];
    LookaheadCell &cell = params.nodes[nodeItem.nodeIndex].lookaheadCell;
    const float count = float(statistics.numSamples);
    cell.radianceFirstMoment += statistics.meanRadiance * count;
    cell.radianceSecondMoment += statistics.meanSquaredRadiance * count;
    cell.radianceWeightedDirectionSum += statistics.meanRadianceWeightedDirection * count;
#if defined(YBI_ENABLE_GUIDING_RGB_STATS)
    cell.rgbRadianceFirstMoment += statistics.meanRadianceRGB * count;
#endif
    cell.numSamples += statistics.numSamples;
}

YBI_GPU float CalculateSplitScore(float meanDelta,
                                  float meanLookahead,
                                  float varianceDelta,
                                  float varianceLookahead)
{
    return (meanDelta - meanLookahead) / Sqrt(varianceDelta + varianceLookahead);
}

YBI_GPU bool ShouldSplitGuidingCell(const PathGuidingKdTreeParams &params,
                                    const LookaheadCell &parent,
                                    const LookaheadCell *lookaheadCells,
                                    uint32_t numLookaheadCells,
                                    float meanRadianceSplitThreshold)
{
    YBI_ASSERT(meanRadianceSplitThreshold >= 0.f);
    if (parent.numSamples < VMM_KD_TREE_MIN_NUM_SAMPLES)
    {
        return false;
    }

    for (uint32_t lookaheadCellIndex = 0; lookaheadCellIndex < numLookaheadCells;
         lookaheadCellIndex++)
    {
        const LookaheadCell &cell = lookaheadCells[lookaheadCellIndex];
        // Section 4.1 requires at least 1000 samples in parent and lookahead.
        // The complement must also be nonempty to define its mean and variance.
        if (cell.numSamples < VMM_KD_TREE_MIN_NUM_SAMPLES || cell.numSamples >= parent.numSamples)
        {
            continue;
        }
        const LookaheadCell delta = CalculateDelta(parent, cell);

        const float oneMinusMeanRadianceSplitThreshold = 1.f - meanRadianceSplitThreshold;
        const float onePlusMeanRadianceSplitThreshold = 1.f + meanRadianceSplitThreshold;

        const float cDeltaPlus = oneMinusMeanRadianceSplitThreshold * delta.numSamples;
        const float cLookaheadPlus =
            parent.numSamples - oneMinusMeanRadianceSplitThreshold * cell.numSamples;

        const float cDeltaMinus = -onePlusMeanRadianceSplitThreshold * delta.numSamples;
        const float cLookaheadMinus =
            onePlusMeanRadianceSplitThreshold * cell.numSamples - parent.numSamples;

        const float meanDelta = CalculateSampleMean(delta);
        const float varianceDelta = CalculateSampleVariance(delta);

        const float meanLookahead = CalculateSampleMean(cell);
        const float varianceLookahead = CalculateSampleVariance(cell);

        const float splitScorePlus = CalculateSplitScore(cDeltaPlus * meanDelta,
                                                         cLookaheadPlus * meanLookahead,
                                                         Sqr(cDeltaPlus) * varianceDelta,
                                                         Sqr(cLookaheadPlus) * varianceLookahead);

        const float splitScoreMinus =
            CalculateSplitScore(cDeltaMinus * meanDelta,
                                cLookaheadMinus * meanLookahead,
                                Sqr(cDeltaMinus) * varianceDelta,
                                Sqr(cLookaheadMinus) * varianceLookahead);

        if (splitScorePlus > VMM_KD_TREE_CRITICAL_Z_SCORE ||
            splitScoreMinus > VMM_KD_TREE_CRITICAL_Z_SCORE)
        {
            return true;
        }

        if (Dot(parent.radianceWeightedDirectionSum, parent.radianceWeightedDirectionSum) == 0.f ||
            Dot(cell.radianceWeightedDirectionSum, cell.radianceWeightedDirectionSum) == 0.f)
            continue;

        const Vec3f parentDirection = CalculateMeanDirection(parent);
        const float parentInvKappa = CalculateMeanDirectionInvKappa(parent);
        const Vec3f lookaheadDirection = CalculateMeanDirection(cell);

        const float invSqrtK = Sqrt(parentInvKappa + CalculateMeanDirectionInvKappa(cell));
        if (invSqrtK >= params.vmfSplitMaxInvSqrtK)
            continue;

        const float position = (invSqrtK / params.vmfSplitMaxInvSqrtK) * (VMFSplitTableSize - 1);
        const float u = (position + 0.5f) / VMFSplitTableSize;
        const float cosBeta = SampleDeviceTexture1D(params.vmfSplitTexture, VMFSplitTableSize, u);
        const float cosine = Clamp(Dot(parentDirection, lookaheadDirection), -1.f, 1.f);
        if (cosine < cosBeta)
            return true;
    }
    return false;
}

YBI_GPU void ShouldSplitGuidingCell(PathGuidingKdTreeParams &params, uint32_t nodeItemIndex)
{
    const PathGuidingKdTreeNode &node = params.nodes[params.nodeItems[nodeItemIndex].nodeIndex];
    if (!node.IsLookahead())
        return;

    const uint32_t ancestorIndex = uint32_t(node.ancestorGuidingCellIndex);
    const LookaheadCell &parent = params.nodes[ancestorIndex].lookaheadCell;
    if (ShouldSplitGuidingCell(
            params, parent, &node.lookaheadCell, 1, params.meanRadianceSplitThreshold))
    {
        const uint32_t bit = 1u << (ancestorIndex & 31u);
        const uint32_t previousMask =
            AtomicOr(&params.guidingCellSplitMask[ancestorIndex >> 5], bit);
        if ((previousMask & bit) == 0)
            AtomicAdd(&params.numGuidingCellSplits, 1u);
    }
}

template <typename Backend>
YBI_GPU void
PromoteGuidingCells(PathGuidingKdTreeParams &params, uint32_t &vmmOffset, uint32_t &frontierOffset)
{
    for (uint32_t nodeIndex = Backend::GetThreadIndex(); nodeIndex < params.numNodes;
         nodeIndex += Backend::GetThreadStride())
    {
        PathGuidingKdTreeNode &node = params.nodes[nodeIndex];
        if (node.IsLookahead() &&
            IsGuidingCellMarkedForSplit(params, uint32_t(node.ancestorGuidingCellIndex)))
            node.lookaheadCell = {};
    }
    Backend::Sync();

    for (uint32_t tileOffset = 0; tileOffset < params.numNodes;
         tileOffset += Backend::GetThreadStride())
    {
        const uint32_t nodeIndex = tileOffset + Backend::GetThreadIndex();
        const bool promote =
            nodeIndex < params.numNodes && IsGuidingCellMarkedForSplit(params, nodeIndex);
        uint32_t childCount = 0;
        if (promote)
        {
            const PathGuidingKdTreeNode &node = params.nodes[nodeIndex];
            YBI_ASSERT(node.IsGuidingCell() && node.HasChildren());
            childCount = uint32_t(params.nodeItems[node.childOffset].sampleCount != 0) +
                         uint32_t(params.nodeItems[node.childOffset + 1].sampleCount != 0);
        }

        const auto vmmScan = Backend::ExclusiveScanAdd(uint32_t(promote));
        const auto frontierScan = Backend::ExclusiveScanAdd(childCount);
        if (Backend::GetThreadIndex() == 0)
        {
            vmmOffset = params.numVMMs;
            params.numVMMs += vmmScan.total;
            frontierOffset = params.numNextNodeItemIndices;
            params.numNextNodeItemIndices += frontierScan.total;
        }
        Backend::Sync();

        if (promote)
        {
            PathGuidingKdTreeNode &node = params.nodes[nodeIndex];
            PathGuidingKdTreeNode &left = params.nodes[node.childOffset];
            PathGuidingKdTreeNode &right = params.nodes[node.childOffset + 1];
            left.guidingCellIndex = node.guidingCellIndex;
            right.guidingCellIndex = int(vmmOffset + vmmScan.prefix);
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
            uint32_t outputIndex = frontierOffset + frontierScan.prefix;
            for (uint32_t childIndex = 0; childIndex < 2; ++childIndex)
            {
                const uint32_t childNodeIndex = uint32_t(node.childOffset) + childIndex;
                if (params.nodeItems[childNodeIndex].sampleCount != 0)
                    params.nextNodeItemIndices[outputIndex++] = childNodeIndex;
            }
        }
        Backend::Sync();
    }

    for (uint32_t nodeIndex = Backend::GetThreadIndex(); nodeIndex < params.numNodes;
         nodeIndex += Backend::GetThreadStride())
    {
        PathGuidingKdTreeNode &node = params.nodes[nodeIndex];
        node.ancestorGuidingCellIndex =
            node.IsGuidingCell() ? -1 : node.GetGuidingAncestorIndex(params.nodes);
    }
}

} // namespace ybi
