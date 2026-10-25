// Implementation adapted from https://github.com/OpenPathGuidingLibrary/openpgl.
// Copyright 2020 Intel Corporation
// SPDX-License-Identifier: Apache-2.0
// Full license: licenses/OpenPGL.txt (relative to the repository root).
#pragma once

#include "core/math.h"
#include "core/octahedral.h"
#include "kernel/guiding/sample.h"
#include "kernel/guiding/types.h"
#include "kernel/guiding/vmf.h"
#include <cfloat>

namespace ybi
{

// Lukas Ruppert, Sebastian Herholz, and Hendrik P. A. Lensch.
// "Robust Fitting of Parallax-Aware Mixtures for Path Guiding", ACM TOG, 2020.
// https://doi.org/10.1145/3386569.3392421

template <typename Backend>
YBI_GPU void UpdateComponentDistances(const PathGuidingKdTreeParams &params,
                                      VMMDistribution &distribution,
                                      SufficientStatistics &statistics,
                                      const VMMUpdateRange &range,
                                      VMMComponentDistanceWorkspace &workspace)
{
    using Float = typename Backend::Float;
    using SampleData = typename Backend::SampleData;
    constexpr uint32_t laneWidth = Backend::SampleWidth;
    const uint32_t numComponents = distribution.numComponents;

    YBI_ASSERT(range.numSamples > 0);
    YBI_ASSERT(numComponents > 0 && numComponents <= VMM_MAX_NUM_COMPONENTS);
    Backend::Sync();
    for (uint32_t k = Backend::GetThreadIndex(); k < numComponents;
         k += Backend::GetThreadStride())
    {
        workspace.batchSumWeights[k] = 0.f;
        workspace.batchSumWeightedInverseDistances[k] = 0.f;
    }
    Backend::Sync();

    // NOTE: calculates the harmonic mean:
    // Σw / Σ(w / d)
    // See Eq. 34 in Ruppert et. al.

    Float batchSumWeights[VMM_MAX_NUM_COMPONENTS]{};
    Float batchSumWeightedInverseDistances[VMM_MAX_NUM_COMPONENTS]{};
    for (uint32_t sampleIndex = Backend::GetThreadIndex() * laneWidth;
         sampleIndex < range.numSamples;
         sampleIndex += Backend::GetThreadStride() * laneWidth)
    {
        const SampleData sample =
            LoadSampleData<SampleData>(params,
                                       range.sampleOffset + sampleIndex,
                                       Min(laneWidth, range.numSamples - sampleIndex));
        YBI_ASSERT(All(sample.pdf > Float(0.f)));
        YBI_ASSERT(All(sample.distance > Float(0.f)));
        const Float sampleWeight = sample.incidentRadianceEstimate / sample.pdf;
        Float assignments[VMM_MAX_NUM_COMPONENTS];
        Float pdf;
        if (!SoftAssignment<Float>(distribution, sample.direction, assignments, pdf))
            continue;
        for (uint32_t k = 0; k < numComponents; ++k)
        {
            if (distribution.weights[k] > FLT_EPSILON)
            {
                const Float componentPdf = assignments[k] * pdf / Float(distribution.weights[k]);
                const Float weight = sampleWeight * assignments[k] * componentPdf;
                batchSumWeights[k] += weight;
                batchSumWeightedInverseDistances[k] += weight / sample.distance;
            }
        }
    }

    for (uint32_t k = 0; k < numComponents; ++k)
    {
        Backend::ReduceAdd(batchSumWeights[k], workspace.batchSumWeights[k]);
        Backend::ReduceAdd(batchSumWeightedInverseDistances[k],
                           workspace.batchSumWeightedInverseDistances[k]);
    }
    Backend::Sync();
    for (uint32_t k = Backend::GetThreadIndex(); k < numComponents;
         k += Backend::GetThreadStride())
    {
        float sumWeightedInverseDistances = workspace.batchSumWeightedInverseDistances[k];
        if (distribution.distance[k] > 0.f)
        {
            sumWeightedInverseDistances +=
                statistics.sumDistanceWeights[k] / distribution.distance[k];
        }
        statistics.sumDistanceWeights[k] += workspace.batchSumWeights[k];
        distribution.distance[k] =
            SafeDivide(statistics.sumDistanceWeights[k], sumWeightedInverseDistances);
    }
    Backend::Sync();
}

template <typename Backend>
YBI_GPU void ReprojectSamples(const PathGuidingKdTreeParams &params, uint32_t updateIndex)
{
    const uint32_t nodeIndex = params.vmmUpdateNodeIndices[updateIndex];
    const PathGuidingKdTreeNode &node = params.nodes[nodeIndex];
    const PathGuidingKdTreeNodeItem &nodeItem = params.nodeItems[nodeIndex];
    YBI_ASSERT(node.IsGuidingCell());
    if (nodeItem.sampleCount == 0)
        return;

    YBI_ASSERT(node.positionStatistics.numSamples > 0);
    const Vec3f pivotPosition = node.positionStatistics.mean;
    const Vec3f variance = node.positionStatistics.variance;
    const float minDistance = 9.f * Sqrt(Sqrt(Max(FLT_EPSILON, Dot(variance, variance))));
    VMMDistribution &distribution = params.vmmDistributions[node.guidingCellIndex];
    if (!distribution.initialized && Backend::GetThreadIndex() == 0)
        distribution.pivotPosition = pivotPosition;

    for (uint32_t sampleOffset = Backend::GetThreadIndex(); sampleOffset < nodeItem.sampleCount;
         sampleOffset += Backend::GetThreadStride())
    {
        const uint32_t outputIndex = nodeItem.sampleOffset + sampleOffset;
        const uint32_t sampleIndex = params.sampleIndices[outputIndex];
        PathGuidingTrainingSample &sample = params.sampleData[sampleIndex];
        Vec3f direction = DecodeOctahedral(sample.octahedralDirection);
        float distance = sample.distance;
        if (distance > 0.f && IsFinite(distance))
        {
            distance = Max(distance, minDistance);
            const Vec3f offset = sample.position - pivotPosition + direction * distance;
            const float shiftedDistance = Length(offset);
            if (shiftedDistance > FLT_EPSILON && IsFinite(shiftedDistance))
            {
                distance = shiftedDistance;
                direction = offset / shiftedDistance;
            }
        }
        sample.octahedralDirection = EncodeOctahedral(direction);
        sample.distance = distance;
    }
}

template <typename Backend>
YBI_GPU void ApplyParallaxShift(const VMMDistribution &distribution,
                                SufficientStatistics &statistics,
                                const Vec3f &shift)
{
    YBI_ASSERT(distribution.numComponents <= VMM_MAX_NUM_COMPONENTS);
    if (Length(shift) < FLT_EPSILON)
        return;

    for (uint32_t k = Backend::GetThreadIndex(); k < distribution.numComponents;
         k += Backend::GetThreadStride())
    {
        const float distance = distribution.distance[k];
        const Vec3f weightedDirection = statistics.sumWeightedDirections[k];
        const float magnitude = Length(weightedDirection);
        const Vec3f offset = SafeDivide(weightedDirection, magnitude) * distance + shift;
        const float shiftedDistance = Length(offset);
        if (distance > 0.f && IsFinite(distance) && magnitude > 0.f && IsFinite(magnitude) &&
            shiftedDistance > FLT_EPSILON && IsFinite(shiftedDistance))
        {
            statistics.sumWeightedDirections[k] = offset * (magnitude / shiftedDistance);
        }
    }
    Backend::Sync();
}

template <typename Backend>
YBI_GPU void PerformRelativeParallaxShift(VMMDistribution &distribution, const Vec3f &shift)
{
    YBI_ASSERT(distribution.numComponents <= VMM_MAX_NUM_COMPONENTS);
    if (Length(shift) < FLT_EPSILON)
        return;

    for (uint32_t k = Backend::GetThreadIndex(); k < distribution.numComponents;
         k += Backend::GetThreadStride())
    {
        const float distance = distribution.distance[k];
        const Vec3f offset = distribution.meanDirection[k] * distance + shift;
        const float shiftedDistance = Length(offset);
        if (distance > 0.f && IsFinite(distance) && IsFinite(shiftedDistance))
        {
            if (shiftedDistance > FLT_EPSILON)
                distribution.meanDirection[k] = offset / shiftedDistance;
            distribution.distance[k] = shiftedDistance;
        }
    }
    Backend::Sync();
    if (Backend::GetThreadIndex() == 0)
        distribution.pivotPosition -= shift;
    Backend::Sync();
}

template <typename Backend>
YBI_GPU void ApplyVMMParallaxShift(const PathGuidingKdTreeParams &params, uint32_t updateIndex)
{
    YBI_ASSERT(updateIndex < params.numVMMUpdateNodes);
    const uint32_t nodeIndex = params.vmmUpdateNodeIndices[updateIndex];
    YBI_ASSERT(nodeIndex < params.numNodes);
    const PathGuidingKdTreeNode &node = params.nodes[nodeIndex];
    YBI_ASSERT(node.IsGuidingCell() && uint32_t(node.guidingCellIndex) < params.numVMMs);
    VMMDistribution &distribution = params.vmmDistributions[node.guidingCellIndex];
    if (!distribution.initialized || node.positionStatistics.numSamples == 0)
        return;

    const Vec3f shift = distribution.pivotPosition - node.positionStatistics.mean;
    ApplyParallaxShift<Backend>(distribution, params.vmmStatistics[node.guidingCellIndex], shift);
    PerformRelativeParallaxShift<Backend>(distribution, shift);
}

} // namespace ybi
