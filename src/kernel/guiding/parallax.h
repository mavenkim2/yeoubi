// Implementation adapted from https://github.com/OpenPathGuidingLibrary/openpgl.
// Copyright 2020 Intel Corporation
// SPDX-License-Identifier: Apache-2.0
// Full license: licenses/OpenPGL.txt (relative to the repository root).
#pragma once

#include "core/math.h"
#include "kernel/guiding/types.h"
#include <cfloat>

namespace ybi
{

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
