#pragma once

#include "core/vec3f.h"
#include "device/cuda_compat.h"
#include "kernel/guiding/sample.h"
#include "kernel/guiding/types.h"
#include "kernel/guiding/vmf.h"
#include "kernel/texture/sampling.h"

namespace ybi
{

#define VMM_KD_TREE_CRITICAL_Z_SCORE 3.7190164855f
#define VMM_KD_TREE_MIN_NUM_SAMPLES 1000u

struct LookaheadCell
{
    float radianceFirstMoment;
    float radianceSecondMoment;
    uint32_t numSamples;
    Vec3f radianceWeightedDirectionSum;

    YBI_GPU Vec3f CalculateMeanDirection() const
    {
        return Normalize(radianceWeightedDirectionSum);
    }

    YBI_GPU float CalculateMeanDirectionInvKappa() const
    {
        const Vec3f r = radianceWeightedDirectionSum / radianceFirstMoment;
        const float meanCosine = Min(Length(r), 1.f);
        if (meanCosine == 1.f)
            return 0.f; // Infinite concentration.

        const float kappa = MeanCosineToKappa(meanCosine);
        const float kappaCothMinusOne = kappa * meanCosine;
        const float yMinusOne =
            (radianceFirstMoment / radianceSecondMoment) * radianceFirstMoment *
            kappaCothMinusOne;
        const float invY = 1.f / (1.f + yMinusOne);
        // Eq. (4), evaluated in 1/y to avoid polynomial overflow.
        const float numerator =
            1.f + invY * (1.69934861f + invY * (5.38753272f + invY * 9.85021305f));
        const float denominator = 1.f + invY * (0.67453491f + invY * 4.31180006f);
        return Sqrt((invY / yMinusOne) * (denominator / numerator));
    }

    YBI_GPU float CalculateSampleMean() const
    {
        return radianceFirstMoment / numSamples;
    }

    YBI_GPU float CalculateSampleVariance() const
    {
        // Variance of the estimated mean, not of an individual sample.
        const double n = numSamples;
        const double b1 = radianceFirstMoment;
        const double numerator = n * radianceSecondMoment - b1 * b1;
        // Accumulated float moments can produce a small negative value through rounding.
        return static_cast<float>(Max(0.0, numerator) / (n * n * n));
    }

    YBI_GPU LookaheadCell CalculateDelta(const LookaheadCell &rhs) const
    {
        YBI_ASSERT(rhs.numSamples <= numSamples);
        LookaheadCell result = {radianceFirstMoment - rhs.radianceFirstMoment,
                                radianceSecondMoment - rhs.radianceSecondMoment,
                                numSamples - rhs.numSamples,
                                radianceWeightedDirectionSum - rhs.radianceWeightedDirectionSum};

        return result;
    }
};

struct GuidingCell
{
};

struct VMMKDTreeWorkspace
{
    float firstMoment;
    float secondMoment;
    Vec3f radianceWeightedDirectionSum;
};

// Clear workspace and synchronize before the first batch; all block threads participate.
template <typename Backend>
YBI_GPU void UpdateLookaheadCell(const GuidingParams &params,
                                 const VMMUpdateRange &range,
                                 VMMKDTreeWorkspace &workspace)
{
    using Float = typename Backend::Float;
    constexpr uint32_t laneWidth = Backend::SampleWidth;

    Float firstMoment(0.f);
    Float secondMoment(0.f);
    Vector3<Float> directionSum{Float(0.f)};

    // TODO: divide by 2 at 8 million samples
    // TODO: samples have to be sorted such that sub ranges are separate lookahead cells
    // TODO: how does openpgl choose split positions again?
    // answer: choose dimension with maximum positional variance... this seems like a pita

    for (uint32_t sampleIndex = Backend::GetThreadIndex() * laneWidth;
         sampleIndex < range.numSamples;
         sampleIndex += Backend::GetThreadStride() * laneWidth)
    {
        const uint32_t numValidSamples = Min(laneWidth, range.numSamples - sampleIndex);
        const Float Li = LoadRadianceContribution<Float>(
            params, range.sampleOffset + sampleIndex, numValidSamples);
        const Vector3<Float> direction = LoadSampleDirection<Float>(
            params, range.sampleOffset + sampleIndex, numValidSamples);
        firstMoment += Li;
        secondMoment += Sqr(Li);
        directionSum += direction * Li;
    }

    Backend::ReduceAdd(firstMoment, workspace.firstMoment);
    Backend::ReduceAdd(secondMoment, workspace.secondMoment);
    Backend::ReduceAdd(directionSum.x, workspace.radianceWeightedDirectionSum.x);
    Backend::ReduceAdd(directionSum.y, workspace.radianceWeightedDirectionSum.y);
    Backend::ReduceAdd(directionSum.z, workspace.radianceWeightedDirectionSum.z);
    Backend::Sync();
}

YBI_GPU float CalculateSplitScore(float meanDelta,
                                  float meanLookahead,
                                  float varianceDelta,
                                  float varianceLookahead)
{
    return (meanDelta - meanLookahead) / Sqrt(varianceDelta + varianceLookahead);
}

YBI_GPU bool SplitGuidingCell(const GuidingParams &params,
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

    const Vec3f parentDirection = parent.CalculateMeanDirection();
    const float parentInvKappa = parent.CalculateMeanDirectionInvKappa();

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
        const LookaheadCell delta = parent.CalculateDelta(cell);

        const float oneMinusMeanRadianceSplitThreshold = 1.f - meanRadianceSplitThreshold;
        const float onePlusMeanRadianceSplitThreshold = 1.f + meanRadianceSplitThreshold;

        const float cDeltaPlus = oneMinusMeanRadianceSplitThreshold * delta.numSamples;
        const float cLookaheadPlus =
            parent.numSamples - oneMinusMeanRadianceSplitThreshold * cell.numSamples;

        const float cDeltaMinus = -onePlusMeanRadianceSplitThreshold * delta.numSamples;
        const float cLookaheadMinus =
            onePlusMeanRadianceSplitThreshold * cell.numSamples - parent.numSamples;

        const float meanDelta = delta.CalculateSampleMean();
        const float varianceDelta = delta.CalculateSampleVariance();

        const float meanLookahead = cell.CalculateSampleMean();
        const float varianceLookahead = cell.CalculateSampleVariance();

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

        const Vec3f lookaheadDirection = cell.CalculateMeanDirection();

        const float invSqrtK = Sqrt(parentInvKappa + cell.CalculateMeanDirectionInvKappa());
        if (invSqrtK >= params.vmfSplitMaxInvSqrtK)
            continue;

        const float position =
            (invSqrtK / params.vmfSplitMaxInvSqrtK) * (VMFSplitTableSize - 1);
        const float u = (position + 0.5f) / VMFSplitTableSize;
        const float cosBeta = SampleDeviceTexture1D(
            params.vmfSplitTexture, VMFSplitTableSize, u);
        const float cosine = Clamp(Dot(parentDirection, lookaheadDirection), -1.f, 1.f);
        if (cosine < cosBeta)
            return true;
    }
    return false;
}

} // namespace ybi
