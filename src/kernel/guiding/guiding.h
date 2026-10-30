// Implementation adapted from https://github.com/OpenPathGuidingLibrary/openpgl.
// Copyright 2020 Intel Corporation
// SPDX-License-Identifier: Apache-2.0
// Full license: licenses/OpenPGL.txt (relative to the repository root).
#pragma once

#include "core/frame.h"
#include "core/vec3f.h"
#include "device/cuda_compat.h"
#include "kernel/guiding/parallax.h"
#include "kernel/guiding/sample.h"
#include "kernel/guiding/types.h"
#include "kernel/guiding/vmf.h"

namespace ybi
{

// Lukas Ruppert, Sebastian Herholz, and Hendrik P. A. Lensch.
// "Robust Fitting of Parallax-Aware Mixtures for Path Guiding", ACM TOG, 2020.
// https://doi.org/10.1145/3386569.3392421

YBI_CPU_GPU void InitializeUniformVMM(VMMDistribution &distribution,
                                      const VMMFitConfiguration &cfg = VMMFitConfiguration{})
{
    YBI_ASSERT(cfg.initialNumComponents > 0 && cfg.initialNumComponents <= VMM_MAX_NUM_COMPONENTS);
    YBI_ASSERT(cfg.initialKappa >= 0.f && IsFinite(cfg.initialKappa));
    distribution = {};
    distribution.split = true;
    distribution.numComponents = cfg.initialNumComponents;
    const float weight = 1.f / float(cfg.initialNumComponents);
    const float goldenRatio = 1.618033988749895f;
    for (uint32_t k = 0; k < cfg.initialNumComponents; ++k)
    {
        const float phi = 2.f * YBI_PI * (float(k) / goldenRatio);
        const float z = 1.f - (2.f * float(k) + 1.f) / float(cfg.initialNumComponents);
        const float radius = Sqrt(Max(0.f, 1.f - z * z));
        distribution.meanDirection[k] = Vec3f(radius * Cos(phi), radius * Sin(phi), z);
        distribution.weights[k] = weight;
        distribution.kappa[k] = cfg.initialKappa;
    }
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
        workspace.currentStats.sumWeights = 0.f;
        workspace.logLikelihood = 0.f;
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
YBI_GPU void NormalizeVMMStatistics(uint32_t numComponents, VMMUpdateWorkspace &workspace)
{
    if (Backend::GetThreadIndex() == 0)
    {
        float sum = 0.f;
        for (uint32_t k = 0; k < numComponents; ++k)
        {
            sum += workspace.currentStats.sumWeightedStats[k];
        }
        workspace.currentStats.sumWeights = sum;
        workspace.normalization =
            sum > FLT_EPSILON ? workspace.currentStats.numSamples / sum : 1.f;
    }
    Backend::Sync();
    for (uint32_t k = Backend::GetThreadIndex(); k < numComponents;
         k += Backend::GetThreadStride())
    {
        workspace.currentStats.sumWeightedDirections[k] *= workspace.normalization;
        workspace.currentStats.sumWeightedStats[k] *= workspace.normalization;
    }
    Backend::Sync();
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
                previousStats.sumDistanceWeights[k] = 0.f;
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
            previousStats.sumWeights += currentStats.sumWeights;
            previousStats.overallNumSamples += currentStats.numSamples;
        }
    }
    Backend::Sync();
}

// All backend participants enter this function together. The caller provides
// shared workspace. The backend supplies worker indices, reductions, sync and reporting.
template <typename Backend, bool Partial = false>
YBI_GPU void UpdateVMM(const PathGuidingKdTreeParams &params,
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
    YBI_ASSERT(range.nonZeroSampleCount <= range.numSamples);
    YBI_ASSERT(range.overallNumSamples >= float(range.nonZeroSampleCount));
    YBI_ASSERT(distribution.numComponents > 0 &&
               distribution.numComponents <= VMM_MAX_NUM_COMPONENTS);
    YBI_ASSERT(cfg.maxEMIterations > 0);
    YBI_ASSERT(cfg.convergenceThreshold >= 0.f);

    if (range.nonZeroSampleCount == 0)
        return;

    float previousLogLikelihood = 0.f;
    for (uint32_t iteration = 0; iteration < cfg.maxEMIterations; ++iteration)
    {
        BeginVMMIteration<Backend>(range.nonZeroSampleCount, distribution.numComponents, workspace);
        VMMStatistics localStats{};
        Float localLogLikelihood(0.f);
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
            const Float sampleWeight = sample.incidentRadianceEstimate / sample.pdf;
            if (sampleWeight == Float(0.f))
                continue;
            Float assignments[VMM_MAX_NUM_COMPONENTS];
            Float pdf;
            if (!SoftAssignment<Float>(distribution, sample.direction, assignments, pdf))
                continue;
            localLogLikelihood += sampleWeight * Log(pdf);
            for (uint32_t k = 0; k < distribution.numComponents; ++k)
            {
                const Float assignedWeight = assignments[k] * sampleWeight;
                localStats.sumWeightedDirections[k] += sample.direction * assignedWeight;
                localStats.sumWeightedStats[k] += assignedWeight;
            }
        }

        ReduceVMMStatistics<Backend>(
            localStats, localLogLikelihood, distribution.numComponents, workspace);
        NormalizeVMMStatistics<Backend>(distribution.numComponents, workspace);
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
    if constexpr (!Partial)
    {
        if (Backend::GetThreadIndex() == 0)
            distribution.initialized = true;
        Backend::Sync();
    }
}

template <typename Backend>
YBI_GPU void CalculateVMMSplitStatistics(const PathGuidingKdTreeParams &params,
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
        workspace.statistics.numSamples[k] = 0.f;
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
        if (sampleWeight == Float(0.f))
            continue;
        Float assignments[VMM_MAX_NUM_COMPONENTS];
        Float V;
        if (!SoftAssignment<Float>(distribution, sample.direction, assignments, V))
            continue;
        const Float directionLengthSquared =
            Sqr(sample.direction.x) + Sqr(sample.direction.y) + Sqr(sample.direction.z);
        Float sumAssignments(0.f);
        for (uint32_t k = 0; k < numComponents; ++k)
            sumAssignments += assignments[k];
        const Float numSamples =
            Select(sumAssignments > Float(0.f),
                   Select(directionLengthSquared > Float(0.f), Float(1.f), Float(0.f)),
                   Float(0.f));
        for (uint32_t k = 0; k < numComponents; ++k)
        {
            const Float assignment = assignments[k];
            localStats.numSamples[k] += numSamples;
            const Float weightedComponentPdf = assignment * V;
            Float chiSquareEstimate = Sqr(Li) * assignment / (Sqr(mcEstimate) * V);
            chiSquareEstimate = chiSquareEstimate - Float(2.f) * assignment * Li / mcEstimate;
            chiSquareEstimate += weightedComponentPdf;
            chiSquareEstimate /= sample.pdf;
            localStats.chiSquareEstimates[k] +=
                Select(directionLengthSquared > Float(0.f), chiSquareEstimate, Float(0.f));

            const Float x = sample.direction.x * Float(workspace.tangents[k].x) +
                            sample.direction.y * Float(workspace.tangents[k].y) +
                            sample.direction.z * Float(workspace.tangents[k].z);
            const Float y = sample.direction.x * Float(workspace.bitangents[k].x) +
                            sample.direction.y * Float(workspace.bitangents[k].y) +
                            sample.direction.z * Float(workspace.bitangents[k].z);
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
        Backend::ReduceAdd(localStats.numSamples[k], workspace.statistics.numSamples[k]);
        Backend::ReduceAdd(localStats.sumWeights[k], workspace.statistics.sumWeights[k]);
        Backend::ReduceAdd(localStats.sumWeightedCovariances[k],
                           workspace.statistics.sumWeightedCovariances[k]);
    }
    Backend::Sync();
    for (uint32_t k = Backend::GetThreadIndex(); k < numComponents;
         k += Backend::GetThreadStride())
    {
        workspace.statistics.chiSquareEstimates[k] =
            Max(0.f,
                SafeDivide(workspace.statistics.chiSquareEstimates[k],
                           workspace.statistics.numSamples[k]));
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
    const float childDistanceWeight = 0.5f * previousStats.sumDistanceWeights[componentIndex];
    const uint32_t childIndex = distribution.numComponents;

    distribution.weights[componentIndex] = childWeight;
    distribution.weights[childIndex] = childWeight;
    distribution.meanDirection[componentIndex] = direction0;
    distribution.meanDirection[childIndex] = direction1;
    distribution.kappa[componentIndex] = childKappa;
    distribution.kappa[childIndex] = childKappa;
    distribution.distance[childIndex] = distribution.distance[componentIndex];
    previousStats.sumWeightedStats[componentIndex] = childStatsWeight;
    previousStats.sumWeightedStats[childIndex] = childStatsWeight;
    previousStats.sumDistanceWeights[componentIndex] = childDistanceWeight;
    previousStats.sumDistanceWeights[childIndex] = childDistanceWeight;
    previousStats.sumWeightedDirections[componentIndex] =
        direction0 * (childMeanCosine * childStatsWeight);
    previousStats.sumWeightedDirections[childIndex] =
        direction1 * (childMeanCosine * childStatsWeight);
    workspace.statistics.chiSquareEstimates[componentIndex] = 0.f;
    workspace.statistics.chiSquareEstimates[childIndex] = 0.f;
    workspace.statistics.numSamples[componentIndex] = 0.f;
    workspace.statistics.numSamples[childIndex] = 0.f;
    workspace.statistics.sumWeights[componentIndex] = 0.f;
    workspace.statistics.sumWeights[childIndex] = 0.f;
    workspace.statistics.sumWeightedCovariances[componentIndex] = Vec3f(0.f);
    workspace.statistics.sumWeightedCovariances[childIndex] = Vec3f(0.f);
    distribution.numComponents++;
    workspace.componentMask |= (1u << componentIndex) | (1u << childIndex);
    return true;
}

template <typename Backend>
YBI_GPU uint32_t SplitVMM(const PathGuidingKdTreeParams &params,
                          VMMDistribution &distribution,
                          SufficientStatistics &previousStats,
                          VMMSplitStatistics<float> &splitStatistics,
                          const VMMUpdateRange &range,
                          float fluxEstimate,
                          const VMMFitConfiguration &cfg,
                          VMMSplitWorkspace &workspace,
                          bool splitComponents = true)
{
    YBI_ASSERT(cfg.splittingThreshold >= 0.f);
    YBI_ASSERT(cfg.map.maxMeanCosine > 0.f);
    CalculateVMMSplitStatistics<Backend>(params, distribution, range, fluxEstimate, workspace);
    for (uint32_t k = Backend::GetThreadIndex(); k < distribution.numComponents;
         k += Backend::GetThreadStride())
    {
        const float batchNumSamples = workspace.statistics.numSamples[k];
        const float numSamples = splitStatistics.numSamples[k] + batchNumSamples;
        workspace.statistics.chiSquareEstimates[k] =
            SafeDivide(splitStatistics.chiSquareEstimates[k] * splitStatistics.numSamples[k] +
                           workspace.statistics.chiSquareEstimates[k] * batchNumSamples,
                       numSamples);
        workspace.statistics.numSamples[k] = numSamples;
        workspace.statistics.sumWeights[k] += splitStatistics.sumWeights[k];
        workspace.statistics.sumWeightedCovariances[k] +=
            splitStatistics.sumWeightedCovariances[k];
    }
    Backend::Sync();
    if (splitComponents && Backend::GetThreadIndex() == 0)
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
    for (uint32_t k = Backend::GetThreadIndex(); k < distribution.numComponents;
         k += Backend::GetThreadStride())
    {
        splitStatistics.chiSquareEstimates[k] = workspace.statistics.chiSquareEstimates[k];
        splitStatistics.numSamples[k] = workspace.statistics.numSamples[k];
        splitStatistics.sumWeights[k] = workspace.statistics.sumWeights[k];
        splitStatistics.sumWeightedCovariances[k] = workspace.statistics.sumWeightedCovariances[k];
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

YBI_CPU_GPU void MergeVMMSplitStatistics(const VMMDistribution &distribution,
                                         VMMSplitStatistics<float> &statistics,
                                         uint32_t component0,
                                         uint32_t component1,
                                         const VMMMergeCandidate &merged)
{
    const float weight0 = distribution.weights[component0];
    const float weight1 = distribution.weights[component1];
    const float sumWeights = statistics.sumWeights[component0] + statistics.sumWeights[component1];
    const Vec3f covariance0 = SafeDivide(statistics.sumWeightedCovariances[component0],
                                         statistics.sumWeights[component0]);
    const Vec3f covariance1 = SafeDivide(statistics.sumWeightedCovariances[component1],
                                         statistics.sumWeights[component1]);
    Vec3f tangent, bitangent;
    BuildOrthonormalBasis(merged.meanDirection, tangent, bitangent);
    const Vec3f direction0 = distribution.meanDirection[component0];
    const Vec3f direction1 = distribution.meanDirection[component1];
    const float angle0 = acosf(Clamp(Dot(direction0, merged.meanDirection), -1.f, 1.f));
    const float angle1 = acosf(Clamp(Dot(direction1, merged.meanDirection), -1.f, 1.f));
    const Vec3f mean0 = Vec3f(Dot(direction0, tangent), Dot(direction0, bitangent), 0.f) *
                        SafeDivide(angle0, Sin(angle0));
    const Vec3f mean1 = Vec3f(Dot(direction1, tangent), Dot(direction1, bitangent), 0.f) *
                        SafeDivide(angle1, Sin(angle1));
    const Vec3f outer0(Sqr(mean0.x), Sqr(mean0.y), mean0.x * mean0.y);
    const Vec3f outer1(Sqr(mean1.x), Sqr(mean1.y), mean1.x * mean1.y);
    statistics.sumWeightedCovariances[component0] =
        SafeDivide((covariance0 + outer0) * weight0 + (covariance1 + outer1) * weight1,
                   merged.weight) *
        sumWeights;
    statistics.sumWeights[component0] = sumWeights;
    statistics.numSamples[component0] = SafeDivide(weight0 * statistics.numSamples[component0] +
                                                       weight1 * statistics.numSamples[component1],
                                                   merged.weight);
    statistics.chiSquareEstimates[component0] += statistics.chiSquareEstimates[component1];

    const uint32_t last = distribution.numComponents - 1;
    statistics.sumWeightedCovariances[component1] = statistics.sumWeightedCovariances[last];
    statistics.sumWeights[component1] = statistics.sumWeights[last];
    statistics.numSamples[component1] = statistics.numSamples[last];
    statistics.chiSquareEstimates[component1] = statistics.chiSquareEstimates[last];
    statistics.sumWeightedCovariances[last] = Vec3f(0.f);
    statistics.sumWeights[last] = 0.f;
    statistics.numSamples[last] = 0.f;
    statistics.chiSquareEstimates[last] = 0.f;
}

YBI_CPU_GPU void MergeVMMComponents(VMMDistribution &distribution,
                                    SufficientStatistics &previousStats,
                                    VMMSplitStatistics<float> &splitStatistics,
                                    uint32_t component0,
                                    uint32_t component1,
                                    const VMMMapConfiguration &cfg,
                                    uint32_t &componentMask)
{
    YBI_ASSERT(component0 < component1 && component1 < distribution.numComponents);
    const VMMMergeCandidate merged =
        CalculateVMMMergedComponent(distribution, component0, component1, cfg);
    MergeVMMSplitStatistics(distribution, splitStatistics, component0, component1, merged);
    const float distance =
        SafeDivide(distribution.weights[component0] * distribution.distance[component0] +
                       distribution.weights[component1] * distribution.distance[component1],
                   merged.weight);
    distribution.weights[component0] = merged.weight;
    distribution.meanDirection[component0] = merged.meanDirection;
    distribution.kappa[component0] = merged.kappa;
    distribution.distance[component0] = distance;
    previousStats.sumWeightedStats[component0] += previousStats.sumWeightedStats[component1];
    previousStats.sumWeightedDirections[component0] +=
        previousStats.sumWeightedDirections[component1];
    previousStats.sumDistanceWeights[component0] += previousStats.sumDistanceWeights[component1];

    const uint32_t last = distribution.numComponents - 1;
    distribution.weights[component1] = distribution.weights[last];
    distribution.meanDirection[component1] = distribution.meanDirection[last];
    distribution.kappa[component1] = distribution.kappa[last];
    distribution.distance[component1] = distribution.distance[last];
    previousStats.sumWeightedStats[component1] = previousStats.sumWeightedStats[last];
    previousStats.sumWeightedDirections[component1] = previousStats.sumWeightedDirections[last];
    previousStats.sumDistanceWeights[component1] = previousStats.sumDistanceWeights[last];
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
    distribution.distance[last] = 0.f;
    previousStats.sumWeightedStats[last] = 0.f;
    previousStats.sumWeightedDirections[last] = Vec3f(0.f);
    previousStats.sumDistanceWeights[last] = 0.f;
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
                          VMMSplitStatistics<float> &splitStatistics,
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
                splitStatistics.numSamples[component0] > 0.f &&
                        splitStatistics.numSamples[component1] > 0.f
                    ? CalculateVMMMergeCost(distribution, component0, component1, cfg.map)
                    : FLOAT_INFINITY;
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
                MergeVMMComponents(distribution,
                                   previousStats,
                                   splitStatistics,
                                   best0,
                                   best1,
                                   cfg.map,
                                   workspace.componentMask);
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

/************************************************************
 * Node dispatch helpers
 ************************************************************/

template <typename Backend>
YBI_GPU void DecayVMMStatistics(VMMDistribution &distribution,
                                SufficientStatistics &statistics,
                                VMMSplitStatistics<float> &splitStatistics,
                                float decay)
{
    for (uint32_t k = Backend::GetThreadIndex(); k < distribution.numComponents;
         k += Backend::GetThreadStride())
    {
        statistics.sumWeightedDirections[k] *= decay;
        statistics.sumWeightedStats[k] *= decay;
        statistics.sumDistanceWeights[k] *= decay;
        splitStatistics.numSamples[k] *= decay;
        splitStatistics.sumWeights[k] *= decay;
        splitStatistics.sumWeightedCovariances[k] *= decay;
    }
    Backend::Sync();
    if (Backend::GetThreadIndex() == 0)
    {
        statistics.numSamples *= decay;
        statistics.sumWeights *= decay;
        distribution.split = false;
    }
    Backend::Sync();
}

template <typename Backend>
YBI_GPU void LoadVMM(const VMMDistribution &storedDistribution,
                     const SufficientStatistics &storedStats,
                     VMMDistribution &distribution,
                     SufficientStatistics &previousStats)
{
    if (Backend::GetThreadIndex() == 0)
    {
        distribution.numComponents = storedDistribution.numComponents;
        distribution.pivotPosition = storedDistribution.pivotPosition;
        distribution.initialized = storedDistribution.initialized;
        distribution.split = storedDistribution.split;
        previousStats.numSamples = storedStats.numSamples;
        previousStats.sumWeights = storedStats.sumWeights;
        previousStats.overallNumSamples = storedStats.overallNumSamples;
        previousStats.numSamplesAfterLastSplit = storedStats.numSamplesAfterLastSplit;
        previousStats.numSamplesAfterLastMerge = storedStats.numSamplesAfterLastMerge;
    }
    Backend::Sync();
    for (uint32_t k = Backend::GetThreadIndex(); k < distribution.numComponents;
         k += Backend::GetThreadStride())
    {
        distribution.weights[k] = storedDistribution.weights[k];
        distribution.meanDirection[k] = storedDistribution.meanDirection[k];
        distribution.kappa[k] = storedDistribution.kappa[k];
        distribution.distance[k] = storedDistribution.distance[k];
        previousStats.sumWeightedDirections[k] = storedStats.sumWeightedDirections[k];
        previousStats.sumWeightedStats[k] = storedStats.sumWeightedStats[k];
        previousStats.sumDistanceWeights[k] = storedStats.sumDistanceWeights[k];
    }
    Backend::Sync();
}

template <typename Backend>
YBI_GPU void StoreVMM(VMMDistribution &storedDistribution,
                      SufficientStatistics &storedStats,
                      const VMMDistribution &distribution,
                      const SufficientStatistics &previousStats)
{
    for (uint32_t k = Backend::GetThreadIndex(); k < distribution.numComponents;
         k += Backend::GetThreadStride())
    {
        storedDistribution.weights[k] = distribution.weights[k];
        storedDistribution.meanDirection[k] = distribution.meanDirection[k];
        storedDistribution.kappa[k] = distribution.kappa[k];
        storedDistribution.distance[k] = distribution.distance[k];
        storedStats.sumWeightedDirections[k] = previousStats.sumWeightedDirections[k];
        storedStats.sumWeightedStats[k] = previousStats.sumWeightedStats[k];
        storedStats.sumDistanceWeights[k] = previousStats.sumDistanceWeights[k];
    }
    if (Backend::GetThreadIndex() == 0)
    {
        storedDistribution.numComponents = distribution.numComponents;
        storedDistribution.pivotPosition = distribution.pivotPosition;
        storedDistribution.initialized = distribution.initialized;
        storedDistribution.split = distribution.split;
        storedStats.numSamples = previousStats.numSamples;
        storedStats.sumWeights = previousStats.sumWeights;
        storedStats.overallNumSamples = previousStats.overallNumSamples;
        storedStats.numSamplesAfterLastSplit = previousStats.numSamplesAfterLastSplit;
        storedStats.numSamplesAfterLastMerge = previousStats.numSamplesAfterLastMerge;
    }
    Backend::Sync();
}

template <typename Backend, bool Partial = false>
YBI_GPU void UpdateNodeVMM(const PathGuidingKdTreeParams &params,
                           uint32_t updateIndex,
                           VMMDistribution &distribution,
                           SufficientStatistics &previousStats,
                           VMMUpdateWorkspace &workspace)
{
    const uint32_t nodeIndex = params.vmmUpdateNodeIndices[updateIndex];
    const PathGuidingKdTreeNode &node = params.nodes[nodeIndex];
    PathGuidingKdTreeNodeItem &nodeItem = params.nodeItems[nodeIndex];
    if constexpr (!Partial)
    {
        const uint32_t nonZeroSampleCount = Backend::template ReduceBlock<uint32_t>(
            [&](uint32_t threadIndex) {
                uint32_t count = 0;
                for (uint32_t sampleOffset = threadIndex; sampleOffset < nodeItem.sampleCount;
                     sampleOffset += VMM_BLOCK_SIZE)
                {
                    const PathGuidingTrainingSample &sample =
                        LoadTrainingSample(params, nodeItem.sampleOffset + sampleOffset);
                    YBI_ASSERT(sample.pdf > 0.f);
                    YBI_ASSERT(sample.GetIncidentRadiance() >= 0.f);
                    count += uint32_t(sample.GetIncidentRadiance() / sample.pdf > 0.f);
                }
                return count;
            });
        if (Backend::GetThreadIndex() == 0)
            nodeItem.nonZeroSampleCount = nonZeroSampleCount;
        Backend::Sync();
    }
    if (nodeItem.nonZeroSampleCount == 0)
        return;
    uint32_t componentMask = ~0u;
    if constexpr (Partial)
    {
        if (!params.vmmFitConfiguration.useSplitAndMerge ||
            !params.vmmFitConfiguration.partialReFit ||
            nodeItem.nonZeroSampleCount < params.vmmFitConfiguration.minSamplesForPartialRefitting)
            return;
        componentMask = params.vmmSplitComponentMasks[node.guidingCellIndex];
        if (componentMask == 0)
            return;
    }
    VMMDistribution &storedDistribution = params.vmmDistributions[node.guidingCellIndex];
    SufficientStatistics &storedStats = params.vmmStatistics[node.guidingCellIndex];
    LoadVMM<Backend>(storedDistribution, storedStats, distribution, previousStats);

    if constexpr (!Partial)
    {
        if (distribution.split)
            DecayVMMStatistics<Backend>(distribution,
                                        previousStats,
                                        params.vmmSplitStatistics[node.guidingCellIndex],
                                        params.vmmFitConfiguration.decayOnSpatialSplit);
        if (distribution.initialized && params.vmmFitConfiguration.useSplitAndMerge &&
            Backend::GetThreadIndex() == 0)
        {
            previousStats.numSamplesAfterLastSplit += nodeItem.nonZeroSampleCount;
            previousStats.numSamplesAfterLastMerge += nodeItem.nonZeroSampleCount;
        }
        Backend::Sync();
    }
    const float overallNumSamples =
        previousStats.overallNumSamples + (Partial ? 0.f : float(nodeItem.nonZeroSampleCount));
    const VMMUpdateRange range{
        nodeItem.sampleOffset, nodeItem.sampleCount, nodeItem.nonZeroSampleCount, overallNumSamples};
    UpdateVMM<Backend, Partial>(params,
                                distribution,
                                previousStats,
                                range,
                                params.vmmFitConfiguration,
                                workspace,
                                componentMask);
    StoreVMM<Backend>(storedDistribution, storedStats, distribution, previousStats);
}

template <typename Backend>
YBI_GPU void SplitNodeVMM(const PathGuidingKdTreeParams &params,
                          uint32_t updateIndex,
                          VMMDistribution &distribution,
                          SufficientStatistics &previousStats,
                          VMMSplitWorkspace &workspace)
{
    const uint32_t nodeIndex = params.vmmUpdateNodeIndices[updateIndex];
    const PathGuidingKdTreeNode &node = params.nodes[nodeIndex];
    const PathGuidingKdTreeNodeItem &nodeItem = params.nodeItems[nodeIndex];
    if (!params.vmmFitConfiguration.useSplitAndMerge || nodeItem.nonZeroSampleCount == 0)
    {
        if (Backend::GetThreadIndex() == 0)
            params.vmmSplitComponentMasks[node.guidingCellIndex] = 0;
        return;
    }
    VMMDistribution &storedDistribution = params.vmmDistributions[node.guidingCellIndex];
    SufficientStatistics &storedStats = params.vmmStatistics[node.guidingCellIndex];
    LoadVMM<Backend>(storedDistribution, storedStats, distribution, previousStats);

    const VMMUpdateRange range{
        nodeItem.sampleOffset, nodeItem.sampleCount, nodeItem.nonZeroSampleCount,
        previousStats.overallNumSamples};
    const bool splitComponents = previousStats.overallNumSamples == float(nodeItem.nonZeroSampleCount) ||
                                 previousStats.numSamplesAfterLastSplit >=
                                     params.vmmFitConfiguration.minSamplesForSplitting;
    const float fluxEstimate = SafeDivide(previousStats.sumWeights, previousStats.numSamples);
    SplitVMM<Backend>(params,
                      distribution,
                      previousStats,
                      params.vmmSplitStatistics[node.guidingCellIndex],
                      range,
                      fluxEstimate,
                      params.vmmFitConfiguration,
                      workspace,
                      splitComponents);
    if (Backend::GetThreadIndex() == 0)
    {
        params.vmmSplitComponentMasks[node.guidingCellIndex] = workspace.componentMask;
        if (splitComponents)
            previousStats.numSamplesAfterLastSplit = 0;
    }
    Backend::Sync();
    StoreVMM<Backend>(storedDistribution, storedStats, distribution, previousStats);
}

template <typename Backend>
YBI_GPU void MergeNodeVMM(const PathGuidingKdTreeParams &params,
                          uint32_t updateIndex,
                          VMMDistribution &distribution,
                          SufficientStatistics &previousStats,
                          VMMMergeWorkspace &workspace)
{
    const uint32_t nodeIndex = params.vmmUpdateNodeIndices[updateIndex];
    const PathGuidingKdTreeNode &node = params.nodes[nodeIndex];
    const PathGuidingKdTreeNodeItem &nodeItem = params.nodeItems[nodeIndex];
    if (!params.vmmFitConfiguration.useSplitAndMerge || nodeItem.nonZeroSampleCount == 0)
        return;
    VMMDistribution &storedDistribution = params.vmmDistributions[node.guidingCellIndex];
    SufficientStatistics &storedStats = params.vmmStatistics[node.guidingCellIndex];
    LoadVMM<Backend>(storedDistribution, storedStats, distribution, previousStats);
    if (previousStats.overallNumSamples != float(nodeItem.nonZeroSampleCount) &&
        previousStats.numSamplesAfterLastMerge < params.vmmFitConfiguration.minSamplesForMerging)
        return;
    MergeVMM<Backend>(distribution,
                      previousStats,
                      params.vmmSplitStatistics[node.guidingCellIndex],
                      params.vmmFitConfiguration,
                      workspace);
    if (Backend::GetThreadIndex() == 0)
        previousStats.numSamplesAfterLastMerge = 0;
    Backend::Sync();
    StoreVMM<Backend>(storedDistribution, storedStats, distribution, previousStats);
}

template <typename Backend>
YBI_GPU void UpdateNodeComponentDistances(const PathGuidingKdTreeParams &params,
                                          uint32_t updateIndex,
                                          VMMDistribution &distribution,
                                          SufficientStatistics &previousStats,
                                          VMMComponentDistanceWorkspace &workspace)
{
    const uint32_t nodeIndex = params.vmmUpdateNodeIndices[updateIndex];
    const PathGuidingKdTreeNode &node = params.nodes[nodeIndex];
    const PathGuidingKdTreeNodeItem &nodeItem = params.nodeItems[nodeIndex];
    if (nodeItem.nonZeroSampleCount == 0)
        return;
    VMMDistribution &storedDistribution = params.vmmDistributions[node.guidingCellIndex];
    SufficientStatistics &storedStats = params.vmmStatistics[node.guidingCellIndex];
    LoadVMM<Backend>(storedDistribution, storedStats, distribution, previousStats);

    const VMMUpdateRange range{
        nodeItem.sampleOffset, nodeItem.sampleCount, nodeItem.nonZeroSampleCount,
        previousStats.overallNumSamples};
    UpdateComponentDistances<Backend>(params, distribution, previousStats, range, workspace);
    StoreVMM<Backend>(storedDistribution, storedStats, distribution, previousStats);
}

} // namespace ybi
