// Implementation adapted from https://github.com/OpenPathGuidingLibrary/openpgl.
// Copyright 2020 Intel Corporation
// SPDX-License-Identifier: Apache-2.0
// Full license: licenses/OpenPGL.txt (relative to the repository root).
#pragma once

#include "core/frame.h"
#include "core/simd_float4.h"
#include "core/vec3f.h"
#include "device/cuda_compat.h"
#include "kernel/guiding/sample.h"
#include "kernel/guiding/types.h"
#include "kernel/guiding/vmf.h"

namespace ybi
{

// Lukas Ruppert, Sebastian Herholz, and Hendrik P. A. Lensch.
// "Robust Fitting of Parallax-Aware Mixtures for Path Guiding", ACM TOG, 2020.
// https://doi.org/10.1145/3386569.3392421

template <typename Float>
YBI_GPU Float SoftAssignment(const VMMDistribution &distribution,
                             const Vector3<Float> &sampleDirection,
                             Float *softAssignment)
{
    Float V(0.f);

    for (uint32_t componentIndex = 0; componentIndex < distribution.numComponents;
         componentIndex++)
    {
        const Float v = Float(distribution.weights[componentIndex]) *
                        EvalVMF<Float>(sampleDirection,
                                       distribution.meanDirection[componentIndex],
                                       distribution.kappa[componentIndex]);

        V += v;

        softAssignment[componentIndex] = v;
    }

    // A zero-density lane has zero assignments and log likelihood log(1) = 0.
    const Float pdf = Select(V == Float(0.f), Float(1.f), V);
    const Float reciprocalPdf = Float(1.f) / pdf;
    for (uint32_t k = 0; k < distribution.numComponents; ++k)
    {
        softAssignment[k] = softAssignment[k] * reciprocalPdf;
    }
    return pdf;
}

template <bool Partial = false>
YBI_CPU_GPU float VMMMapComponentWeight(uint32_t numComponents,
                                        uint32_t componentIndex,
                                        const SufficientStatistics &currentStats,
                                        const SufficientStatistics &previousStats,
                                        const VMMMapConfiguration &cfg)
{
    const float previousNumSamples = Partial ? 0.f : previousStats.numSamples;
    const float previousWeight = Partial ? 0.f : previousStats.sumWeightedStats[componentIndex];
    const float numSamples = currentStats.numSamples + previousNumSamples;
    return SafeDivide(cfg.weightPrior + currentStats.sumWeightedStats[componentIndex] +
                          previousWeight,
                      cfg.weightPrior * numComponents + numSamples);
}

template <bool Partial = false>
YBI_CPU_GPU void VMMMapStateUpdate(VMMDistribution &distribution,
                                   uint32_t componentIndex,
                                   const SufficientStatistics &currentStats,
                                   const SufficientStatistics &previousStats,
                                   float overallNumSamples,
                                   const VMMMapConfiguration &cfg = VMMMapConfiguration{},
                                   float weightScale = 1.f)
{
    YBI_ASSERT(distribution.numComponents <= VMM_MAX_NUM_COMPONENTS);
    if (componentIndex >= distribution.numComponents)
    {
        return;
    }

    const float previousNumSamples = Partial ? 0.f : previousStats.numSamples;
    const float numSamples = currentStats.numSamples + previousNumSamples;
    YBI_ASSERT(numSamples > 0.f);

    const float newWeight =
        VMMMapComponentWeight<Partial>(
            distribution.numComponents, componentIndex, currentStats, previousStats, cfg) *
        weightScale;

    const Vec3f currentMean = SafeDivide(currentStats.sumWeightedDirections[componentIndex],
                                         currentStats.sumWeightedStats[componentIndex]);
    const Vec3f previousMean =
        Partial ? Vec3f(0.f)
                : SafeDivide(previousStats.sumWeightedDirections[componentIndex],
                             previousStats.sumWeightedStats[componentIndex]);
    const float currentEstimationWeight = currentStats.numSamples / numSamples;
    const Vec3f mean =
        currentMean * currentEstimationWeight + previousMean * (1.f - currentEstimationWeight);
    const float meanCosine = Length(mean);
    const Vec3f newDirection =
        SafeDivide(mean, meanCosine, distribution.meanDirection[componentIndex]);

    const float partialNumSamples = newWeight * overallNumSamples;
    const float regularizedMeanCosine = SafeDivide(
        cfg.meanCosinePrior * cfg.meanCosinePriorStrength + meanCosine * partialNumSamples,
        cfg.meanCosinePriorStrength + partialNumSamples,
        meanCosine);
    const float newMeanCosine =
        Clamp(regularizedMeanCosine, 0.f, Min(cfg.maxMeanCosine, FLOAT_BEFORE_ONE));

    distribution.weights[componentIndex] = newWeight;
    distribution.meanDirection[componentIndex] = newDirection;
    distribution.kappa[componentIndex] = MeanCosineToKappa(newMeanCosine);
}

template <typename Backend>
YBI_GPU void
BeginVMMIteration(uint32_t numSamples, uint32_t numComponents, VMMUpdateWorkspace &workspace)
{
    Backend::Sync();
    for (uint32_t k = Backend::GetThreadIndex(); k < numComponents;
         k += Backend::GetThreadStride())
    {
        workspace.currentStats.sumWeightedDirections[k] = Vec3f(0.f);
        workspace.currentStats.sumWeightedStats[k] = 0.f;
    }
    if (Backend::GetThreadIndex() == 0)
    {
        workspace.currentStats.numSamples = float(numSamples);
        workspace.logLikelihood = 0.f;
        workspace.sumWeights = 0.f;
        workspace.normalization = 0.f;
    }
    Backend::Sync();
}

template <typename Backend>
YBI_GPU void ReduceVMMStatistics(const typename Backend::VMMStatistics &localStats,
                                 typename Backend::Float localLogLikelihood,
                                 uint32_t numComponents,
                                 VMMUpdateWorkspace &workspace)
{
    for (uint32_t k = 0; k < numComponents; ++k)
    {
        Backend::ReduceAdd(localStats.sumWeightedDirections[k],
                           workspace.currentStats.sumWeightedDirections[k]);
        Backend::ReduceAdd(localStats.sumWeightedStats[k],
                           workspace.currentStats.sumWeightedStats[k]);
    }
    Backend::ReduceAdd(localLogLikelihood, workspace.logLikelihood);
    Backend::Sync();
}

template <typename Backend>
YBI_GPU bool NormalizeVMMStatistics(uint32_t numComponents, VMMUpdateWorkspace &workspace)
{
    if (Backend::GetThreadIndex() == 0)
    {
        float sum = 0.f;
        for (uint32_t k = 0; k < numComponents; ++k)
        {
            sum += workspace.currentStats.sumWeightedStats[k];
        }
        workspace.sumWeights = sum;
        workspace.normalization = SafeDivide(workspace.currentStats.numSamples, sum);
    }
    Backend::Sync();
    if (workspace.sumWeights == 0.f)
    {
        return false;
    }
    for (uint32_t k = Backend::GetThreadIndex(); k < numComponents;
         k += Backend::GetThreadStride())
    {
        workspace.currentStats.sumWeightedDirections[k] *= workspace.normalization;
        workspace.currentStats.sumWeightedStats[k] *= workspace.normalization;
    }
    Backend::Sync();
    return true;
}

template <typename Backend, bool Partial = false>
YBI_GPU void UpdateVMMComponents(VMMDistribution &distribution,
                                 VMMUpdateWorkspace &workspace,
                                 const SufficientStatistics &previousStats,
                                 float overallNumSamples,
                                 const VMMMapConfiguration &cfg,
                                 uint32_t componentMask = ~0u)
{
    const SufficientStatistics &currentStats = workspace.currentStats;
    float weightScale = 1.f;
    if constexpr (Partial)
    {
        if (Backend::GetThreadIndex() == 0)
        {
            float fixedWeight = 0.f;
            float selectedWeight = 0.f;
            for (uint32_t k = 0; k < distribution.numComponents; ++k)
            {
                if (componentMask & (1u << k))
                {
                    selectedWeight += VMMMapComponentWeight<true>(
                        distribution.numComponents, k, currentStats, previousStats, cfg);
                }
                else
                {
                    fixedWeight += distribution.weights[k];
                }
            }
            workspace.partialWeightScale = SafeDivide(Max(0.f, 1.f - fixedWeight), selectedWeight);
        }
        Backend::Sync();
        weightScale = workspace.partialWeightScale;
    }
    for (uint32_t componentIndex = Backend::GetThreadIndex();
         componentIndex < distribution.numComponents;
         componentIndex += Backend::GetThreadStride())
    {
        if constexpr (Partial)
        {
            if ((componentMask & (1u << componentIndex)) == 0 || weightScale == 0.f)
            {
                continue;
            }
        }
        VMMMapStateUpdate<Partial>(distribution,
                                   componentIndex,
                                   currentStats,
                                   previousStats,
                                   overallNumSamples,
                                   cfg,
                                   weightScale);
    }
    Backend::Sync();
}

template <typename Backend, bool Partial = false>
YBI_GPU void CommitVMMStatistics(SufficientStatistics &previousStats,
                                 const SufficientStatistics &currentStats,
                                 uint32_t numComponents,
                                 uint32_t componentMask = ~0u)
{
    for (uint32_t k = Backend::GetThreadIndex(); k < numComponents;
         k += Backend::GetThreadStride())
    {
        if constexpr (Partial)
        {
            if (componentMask & (1u << k))
            {
                previousStats.sumWeightedDirections[k] = currentStats.sumWeightedDirections[k];
                previousStats.sumWeightedStats[k] = currentStats.sumWeightedStats[k];
            }
        }
        else
        {
            previousStats.sumWeightedDirections[k] += currentStats.sumWeightedDirections[k];
            previousStats.sumWeightedStats[k] += currentStats.sumWeightedStats[k];
        }
    }
    if constexpr (Partial)
    {
        Backend::Sync();
    }
    if (Backend::GetThreadIndex() == 0)
    {
        if constexpr (Partial)
        {
            float numSamples = 0.f;
            for (uint32_t k = 0; k < numComponents; ++k)
            {
                numSamples += previousStats.sumWeightedStats[k];
            }
            previousStats.numSamples = numSamples;
        }
        else
        {
            previousStats.numSamples += currentStats.numSamples;
        }
    }
    Backend::Sync();
}

// All backend participants enter this function together. The caller provides
// shared workspace. The backend supplies worker indices, reductions, sync and reporting.
template <typename Backend, bool Partial = false>
YBI_GPU void UpdateVMM(const GuidingParams &params,
                       VMMDistribution &distribution,
                       SufficientStatistics &previousStats,
                       const VMMUpdateRange &range,
                       const VMMFitConfiguration &cfg,
                       VMMUpdateWorkspace &workspace,
                       uint32_t componentMask = ~0u)
{
    using Float = typename Backend::Float;
    using SampleData = typename Backend::SampleData;
    using VMMStatistics = typename Backend::VMMStatistics;
    constexpr uint32_t laneWidth = Backend::SampleWidth;

    YBI_ASSERT(range.numSamples > 0);
    YBI_ASSERT(range.overallNumSamples >= float(range.numSamples));
    YBI_ASSERT(distribution.numComponents > 0 &&
               distribution.numComponents <= VMM_MAX_NUM_COMPONENTS);
    YBI_ASSERT(cfg.maxEMIterations > 0);
    YBI_ASSERT(cfg.convergenceThreshold >= 0.f);

    float previousLogLikelihood = 0.f;
    for (uint32_t iteration = 0; iteration < cfg.maxEMIterations; ++iteration)
    {
        BeginVMMIteration<Backend>(range.numSamples, distribution.numComponents, workspace);
        VMMStatistics localStats{};
        Float localLogLikelihood(0.f);
        for (uint32_t sampleIndex = Backend::GetThreadIndex() * laneWidth;
             sampleIndex < range.numSamples;
             sampleIndex += Backend::GetThreadStride() * laneWidth)
        {
            const SampleData sample = LoadSampleData<SampleData>(
                params, range.sampleOffset + sampleIndex, range.numSamples);
            YBI_ASSERT(All(sample.pdf > Float(0.f)));
            YBI_ASSERT(All(sample.incidentRadianceEstimate >= Float(0.f)));
            const Float sampleWeight = sample.incidentRadianceEstimate / sample.pdf;
            Float assignments[VMM_MAX_NUM_COMPONENTS];
            const Float pdf =
                SoftAssignment<Float>(distribution, sample.sampleDirection, assignments);
            localLogLikelihood += sampleWeight * Log(pdf);
            for (uint32_t k = 0; k < distribution.numComponents; ++k)
            {
                const Float assignedWeight = assignments[k] * sampleWeight;
                localStats.sumWeightedDirections[k] += sample.sampleDirection * assignedWeight;
                localStats.sumWeightedStats[k] += assignedWeight;
            }
        }

        ReduceVMMStatistics<Backend>(
            localStats, localLogLikelihood, distribution.numComponents, workspace);
        if (!NormalizeVMMStatistics<Backend>(distribution.numComponents, workspace))
        {
            return;
        }
        UpdateVMMComponents<Backend, Partial>(distribution,
                                              workspace,
                                              previousStats,
                                              range.overallNumSamples,
                                              cfg.map,
                                              componentMask);

        const float logLikelihood = workspace.logLikelihood;
        Backend::ReportIteration(
            iteration, logLikelihood, previousLogLikelihood, cfg.printLogLikelihood);
        if (iteration > 0 && Abs(logLikelihood - previousLogLikelihood) <=
                                 cfg.convergenceThreshold * Abs(previousLogLikelihood))
        {
            break;
        }
        previousLogLikelihood = logLikelihood;
    }
    // Commit only the final E-step statistics, not every EM iteration.
    CommitVMMStatistics<Backend, Partial>(
        previousStats, workspace.currentStats, distribution.numComponents, componentMask);
}

template <typename Backend>
YBI_GPU void CalculateVMMSplitStatistics(const GuidingParams &params,
                                         const VMMDistribution &distribution,
                                         const VMMUpdateRange &range,
                                         float fluxEstimate,
                                         VMMSplitWorkspace &workspace)
{
    using Float = typename Backend::Float;
    using SampleData = typename Backend::SampleData;
    constexpr uint32_t laneWidth = Backend::SampleWidth;

    YBI_ASSERT(range.numSamples > 0);
    YBI_ASSERT(fluxEstimate >= 0.f);
    YBI_ASSERT(distribution.numComponents > 0 &&
               distribution.numComponents <= VMM_MAX_NUM_COMPONENTS);
    Backend::Sync();
    const uint32_t numComponents = distribution.numComponents;
    for (uint32_t k = Backend::GetThreadIndex(); k < numComponents;
         k += Backend::GetThreadStride())
    {
        workspace.statistics.chiSquareEstimates[k] = 0.f;
        workspace.statistics.sumWeights[k] = 0.f;
        workspace.statistics.sumWeightedCovariances[k] = Vec3f(0.f);
        BuildOrthonormalBasis(
            distribution.meanDirection[k], workspace.tangents[k], workspace.bitangents[k]);
    }
    if (Backend::GetThreadIndex() == 0)
    {
        workspace.numSplits = 0;
        workspace.componentMask = 0;
    }
    Backend::Sync();

    if (fluxEstimate == 0.f)
    {
        return;
    }

    VMMSplitStatistics<Float> localStats{};
    const Float mcEstimate(fluxEstimate);
    for (uint32_t sampleIndex = Backend::GetThreadIndex() * laneWidth;
         sampleIndex < range.numSamples;
         sampleIndex += Backend::GetThreadStride() * laneWidth)
    {
        const SampleData sample =
            LoadSampleData<SampleData>(params,
                                       range.sampleOffset + sampleIndex,
                                       Min(laneWidth, range.numSamples - sampleIndex));
        YBI_ASSERT(All(sample.pdf > Float(0.f)));
        YBI_ASSERT(All(sample.incidentRadianceEstimate >= Float(0.f)));
        const Float Li = sample.incidentRadianceEstimate;
        const Float sampleWeight = Li / sample.pdf;
        Float assignments[VMM_MAX_NUM_COMPONENTS];
        const Float V = SoftAssignment<Float>(distribution, sample.sampleDirection, assignments);
        const Float directionLengthSquared = Sqr(sample.sampleDirection.x) +
                                             Sqr(sample.sampleDirection.y) +
                                             Sqr(sample.sampleDirection.z);
        for (uint32_t k = 0; k < numComponents; ++k)
        {
            const Float assignment = assignments[k];
            const Float weightedComponentPdf = assignment * V;
            Float chiSquareEstimate = Sqr(Li) * assignment / (Sqr(mcEstimate) * V);
            chiSquareEstimate = chiSquareEstimate - Float(2.f) * assignment * Li / mcEstimate;
            chiSquareEstimate += weightedComponentPdf;
            chiSquareEstimate /= sample.pdf;
            localStats.chiSquareEstimates[k] +=
                Select(directionLengthSquared > Float(0.f), chiSquareEstimate, Float(0.f));

            const Float x = sample.sampleDirection.x * Float(workspace.tangents[k].x) +
                            sample.sampleDirection.y * Float(workspace.tangents[k].y) +
                            sample.sampleDirection.z * Float(workspace.tangents[k].z);
            const Float y = sample.sampleDirection.x * Float(workspace.bitangents[k].x) +
                            sample.sampleDirection.y * Float(workspace.bitangents[k].y) +
                            sample.sampleDirection.z * Float(workspace.bitangents[k].z);
            const Float assignedWeight = assignment * sampleWeight;
            localStats.sumWeights[k] += assignedWeight;
            localStats.sumWeightedCovariances[k] +=
                Vector3<Float>(Sqr(x), Sqr(y), x * y) * assignedWeight;
        }
    }
    for (uint32_t k = 0; k < numComponents; ++k)
    {
        Backend::ReduceAdd(localStats.chiSquareEstimates[k],
                           workspace.statistics.chiSquareEstimates[k]);
        Backend::ReduceAdd(localStats.sumWeights[k], workspace.statistics.sumWeights[k]);
        Backend::ReduceAdd(localStats.sumWeightedCovariances[k],
                           workspace.statistics.sumWeightedCovariances[k]);
    }
    Backend::Sync();
    for (uint32_t k = Backend::GetThreadIndex(); k < numComponents;
         k += Backend::GetThreadStride())
    {
        workspace.statistics.chiSquareEstimates[k] =
            Max(0.f, workspace.statistics.chiSquareEstimates[k] / float(range.numSamples));
    }
    Backend::Sync();
}

YBI_CPU_GPU bool SplitVMMComponent(VMMDistribution &distribution,
                                   SufficientStatistics &previousStats,
                                   uint32_t componentIndex,
                                   VMMSplitWorkspace &workspace,
                                   const VMMFitConfiguration &cfg)
{
    YBI_ASSERT(componentIndex < distribution.numComponents);
    if (distribution.numComponents == VMM_MAX_NUM_COMPONENTS ||
        workspace.statistics.sumWeights[componentIndex] == 0.f)
    {
        return false;
    }
    const float maxMeanCosine = Min(cfg.map.maxMeanCosine, FLOAT_BEFORE_ONE);
    if (distribution.kappa[componentIndex] >= 0.9f * MeanCosineToKappa(maxMeanCosine))
    {
        return false;
    }
    const Vec3f covariance = workspace.statistics.sumWeightedCovariances[componentIndex] /
                             workspace.statistics.sumWeights[componentIndex];
    const float D = 0.5f * Sqrt(Sqr(covariance.x - covariance.y) + 4.f * Sqr(covariance.z));
    if (D <= 1e-8f)
    {
        return false;
    }
    const float eigenValue = 0.5f * (covariance.x + covariance.y) + D;
    Vec3f eigenVector(covariance.z, eigenValue - covariance.x, 0.f);
    if (LengthSquared(eigenVector) == 0.f)
    {
        eigenVector = Vec3f(1.f, 0.f, 0.f);
    }
    eigenVector = Normalize(eigenVector);

    const float angle = 0.5f * eigenValue;
    const float cosine = Cos(angle);
    const Vec3f meanDirection = distribution.meanDirection[componentIndex];
    const Vec3f tangentDirection = workspace.tangents[componentIndex] * eigenVector.x +
                                   workspace.bitangents[componentIndex] * eigenVector.y;
    const Vec3f offset = tangentDirection * Sin(angle);
    const Vec3f direction0 = Normalize(meanDirection * cosine + offset);
    const Vec3f direction1 = Normalize(meanDirection * cosine - offset);
    const float meanCosine = KappaToMeanCosine(distribution.kappa[componentIndex]);
    const float childMeanCosine =
        Min(maxMeanCosine, meanCosine / Abs(Dot(meanDirection, direction0)));
    const float childKappa = MeanCosineToKappa(childMeanCosine);
    const float childWeight = 0.5f * distribution.weights[componentIndex];
    const float childStatsWeight = 0.5f * previousStats.sumWeightedStats[componentIndex];
    const uint32_t childIndex = distribution.numComponents;

    distribution.weights[componentIndex] = childWeight;
    distribution.weights[childIndex] = childWeight;
    distribution.meanDirection[componentIndex] = direction0;
    distribution.meanDirection[childIndex] = direction1;
    distribution.kappa[componentIndex] = childKappa;
    distribution.kappa[childIndex] = childKappa;
    previousStats.sumWeightedStats[componentIndex] = childStatsWeight;
    previousStats.sumWeightedStats[childIndex] = childStatsWeight;
    previousStats.sumWeightedDirections[componentIndex] =
        direction0 * (childMeanCosine * childStatsWeight);
    previousStats.sumWeightedDirections[childIndex] =
        direction1 * (childMeanCosine * childStatsWeight);
    workspace.statistics.chiSquareEstimates[componentIndex] = 0.f;
    workspace.statistics.chiSquareEstimates[childIndex] = 0.f;
    workspace.statistics.sumWeights[componentIndex] = 0.f;
    workspace.statistics.sumWeights[childIndex] = 0.f;
    workspace.statistics.sumWeightedCovariances[componentIndex] = Vec3f(0.f);
    workspace.statistics.sumWeightedCovariances[childIndex] = Vec3f(0.f);
    distribution.numComponents++;
    workspace.componentMask |= (1u << componentIndex) | (1u << childIndex);
    return true;
}

template <typename Backend>
YBI_GPU uint32_t SplitVMM(const GuidingParams &params,
                          VMMDistribution &distribution,
                          SufficientStatistics &previousStats,
                          const VMMUpdateRange &range,
                          float fluxEstimate,
                          const VMMFitConfiguration &cfg,
                          VMMSplitWorkspace &workspace)
{
    YBI_ASSERT(cfg.splittingThreshold >= 0.f);
    YBI_ASSERT(cfg.map.maxMeanCosine > 0.f);
    CalculateVMMSplitStatistics<Backend>(params, distribution, range, fluxEstimate, workspace);
    if (Backend::GetThreadIndex() == 0)
    {
        const uint32_t numComponents = distribution.numComponents;
        uint32_t candidates[VMM_MAX_NUM_COMPONENTS];
        for (uint32_t k = 0; k < numComponents; ++k)
        {
            uint32_t position = k;
            while (position > 0 &&
                   workspace.statistics.chiSquareEstimates[candidates[position - 1]] <
                       workspace.statistics.chiSquareEstimates[k])
            {
                candidates[position] = candidates[position - 1];
                --position;
            }
            candidates[position] = k;
        }
        for (uint32_t i = 0;
             i < numComponents && distribution.numComponents < VMM_MAX_NUM_COMPONENTS;
             ++i)
        {
            const uint32_t k = candidates[i];
            if (workspace.statistics.chiSquareEstimates[k] > cfg.splittingThreshold &&
                SplitVMMComponent(distribution, previousStats, k, workspace, cfg))
            {
                ++workspace.numSplits;
            }
        }
    }
    Backend::Sync();
    return workspace.numSplits;
}

YBI_CPU_GPU VMMMergeCandidate CalculateVMMMergedComponent(const VMMDistribution &distribution,
                                                          uint32_t component0,
                                                          uint32_t component1,
                                                          const VMMMapConfiguration &cfg)
{
    const float weight0 = distribution.weights[component0];
    const float weight1 = distribution.weights[component1];
    VMMMergeCandidate merged;
    merged.weight = weight0 + weight1;
    const Vec3f moment = distribution.meanDirection[component0] *
                             (weight0 * KappaToMeanCosine(distribution.kappa[component0])) +
                         distribution.meanDirection[component1] *
                             (weight1 * KappaToMeanCosine(distribution.kappa[component1]));
    const float momentLength = Length(moment);
    merged.meanDirection =
        momentLength > 0.f ? moment / momentLength : distribution.meanDirection[component0];
    const float meanCosine =
        Min(SafeDivide(momentLength, merged.weight), Min(cfg.maxMeanCosine, FLOAT_BEFORE_ONE));
    merged.kappa = MeanCosineToKappa(meanCosine);
    if (merged.kappa < 1e-3f)
    {
        merged.kappa = 0.f;
    }
    return merged;
}

YBI_CPU_GPU float IntegratedVMFProductRatio(const Vec3f &direction0,
                                            float kappa0,
                                            const Vec3f &direction1,
                                            float kappa1,
                                            const VMMMergeCandidate &merged)
{
    const double x = double(direction0.x) * kappa0 + double(direction1.x) * kappa1 -
                     double(merged.meanDirection.x) * merged.kappa;
    const double y = double(direction0.y) * kappa0 + double(direction1.y) * kappa1 -
                     double(merged.meanDirection.y) * merged.kappa;
    const double z = double(direction0.z) * kappa0 + double(direction1.z) * kappa1 -
                     double(merged.meanDirection.z) * merged.kappa;
    const double kappa = sqrt(x * x + y * y + z * z);
    return (VMFNormalization(kappa0) / VMFNormalization(merged.kappa)) *
           (VMFNormalization(kappa1) / VMFNormalization(float(kappa))) *
           Exp(float(kappa - kappa0 - kappa1 + merged.kappa));
}

YBI_CPU_GPU float CalculateVMMMergeCost(const VMMDistribution &distribution,
                                        uint32_t component0,
                                        uint32_t component1,
                                        const VMMMapConfiguration &cfg)
{
    const VMMMergeCandidate merged =
        CalculateVMMMergedComponent(distribution, component0, component1, cfg);
    if (merged.weight == 0.f)
    {
        return 0.f;
    }
    const float weight0 = distribution.weights[component0];
    const float weight1 = distribution.weights[component1];
    const Vec3f direction0 = distribution.meanDirection[component0];
    const Vec3f direction1 = distribution.meanDirection[component1];
    const float kappa0 = distribution.kappa[component0];
    const float kappa1 = distribution.kappa[component1];
    float integral = 0.f;
    if (weight0 > 0.f)
    {
        integral += Sqr(weight0) *
                    IntegratedVMFProductRatio(direction0, kappa0, direction0, kappa0, merged);
    }
    if (weight1 > 0.f)
    {
        integral += Sqr(weight1) *
                    IntegratedVMFProductRatio(direction1, kappa1, direction1, kappa1, merged);
    }
    if (weight0 > 0.f && weight1 > 0.f)
    {
        integral += 2.f * weight0 * weight1 *
                    IntegratedVMFProductRatio(direction0, kappa0, direction1, kappa1, merged);
    }
    return Max(0.f, integral / merged.weight - merged.weight);
}

YBI_CPU_GPU void MergeVMMComponents(VMMDistribution &distribution,
                                    SufficientStatistics &previousStats,
                                    uint32_t component0,
                                    uint32_t component1,
                                    const VMMMapConfiguration &cfg,
                                    uint32_t &componentMask)
{
    YBI_ASSERT(component0 < component1 && component1 < distribution.numComponents);
    const VMMMergeCandidate merged =
        CalculateVMMMergedComponent(distribution, component0, component1, cfg);
    distribution.weights[component0] = merged.weight;
    distribution.meanDirection[component0] = merged.meanDirection;
    distribution.kappa[component0] = merged.kappa;
    previousStats.sumWeightedStats[component0] += previousStats.sumWeightedStats[component1];
    previousStats.sumWeightedDirections[component0] +=
        previousStats.sumWeightedDirections[component1];

    const uint32_t last = distribution.numComponents - 1;
    distribution.weights[component1] = distribution.weights[last];
    distribution.meanDirection[component1] = distribution.meanDirection[last];
    distribution.kappa[component1] = distribution.kappa[last];
    previousStats.sumWeightedStats[component1] = previousStats.sumWeightedStats[last];
    previousStats.sumWeightedDirections[component1] = previousStats.sumWeightedDirections[last];
    const bool movedComponentSelected = (componentMask & (1u << last)) != 0;
    componentMask &= ~((1u << component1) | (1u << last));
    if (component1 != last && movedComponentSelected)
    {
        componentMask |= 1u << component1;
    }
    componentMask |= 1u << component0;

    distribution.weights[last] = 0.f;
    distribution.meanDirection[last] = Vec3f(0.f);
    distribution.kappa[last] = 0.f;
    previousStats.sumWeightedStats[last] = 0.f;
    previousStats.sumWeightedDirections[last] = Vec3f(0.f);
    --distribution.numComponents;
}

YBI_CPU_GPU void VMMComponentPairFromIndex(uint32_t pair,
                                           uint32_t numComponents,
                                           uint32_t &component0,
                                           uint32_t &component1)
{
    const float d = 2.f * numComponents - 1.f;
    component0 = uint32_t(0.5f * (d - Sqrt(d * d - 8.f * pair)));
    const uint32_t rowStart = component0 * (2 * numComponents - component0 - 1) / 2;
    component1 = component0 + 1 + pair - rowStart;
}

template <typename Backend>
YBI_GPU uint32_t MergeVMM(VMMDistribution &distribution,
                          SufficientStatistics &previousStats,
                          const VMMFitConfiguration &cfg,
                          VMMMergeWorkspace &workspace)
{
    Backend::Sync();
    YBI_ASSERT(distribution.numComponents > 0 &&
               distribution.numComponents <= VMM_MAX_NUM_COMPONENTS);
    YBI_ASSERT(cfg.mergingThreshold >= 0.f);
    YBI_ASSERT(cfg.map.maxMeanCosine > 0.f);
    if (Backend::GetThreadIndex() == 0)
    {
        workspace.numMerges = 0;
        workspace.componentMask = 0;
    }
    Backend::Sync();
    while (distribution.numComponents > 1)
    {
        const uint32_t numComponents = distribution.numComponents;
        const uint32_t numPairs = numComponents * (numComponents - 1) / 2;
        for (uint32_t pair = Backend::GetThreadIndex(); pair < numPairs;
             pair += Backend::GetThreadStride())
        {
            uint32_t component0, component1;
            VMMComponentPairFromIndex(pair, numComponents, component0, component1);
            workspace.costs[pair] =
                CalculateVMMMergeCost(distribution, component0, component1, cfg.map);
        }
        Backend::Sync();
        if (Backend::GetThreadIndex() == 0)
        {
            float minCost = cfg.mergingThreshold;
            uint32_t best0 = 0;
            uint32_t best1 = 0;
            uint32_t pair = 0;
            for (uint32_t component0 = 0; component0 < numComponents; ++component0)
            {
                for (uint32_t component1 = component0 + 1; component1 < numComponents;
                     ++component1, ++pair)
                {
                    if (workspace.costs[pair] < minCost)
                    {
                        minCost = workspace.costs[pair];
                        best0 = component0;
                        best1 = component1;
                    }
                }
            }
            workspace.merged = best0 != best1;
            if (workspace.merged)
            {
                MergeVMMComponents(
                    distribution, previousStats, best0, best1, cfg.map, workspace.componentMask);
                ++workspace.numMerges;
            }
        }
        Backend::Sync();
        if (!workspace.merged)
        {
            break;
        }
    }
    return workspace.numMerges;
}

} // namespace ybi
