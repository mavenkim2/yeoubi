#pragma once

#include "core/vec3f.h"
#include "device/cuda_compat.h"
#include "kernel/guiding/vmf.h"

namespace ybi
{

// 1. calculate fluence estimate, which is just 1/N * sum of mcEstimate
// 2. calculate softassignment weight per sample per vmm component. multiply by sample weight
// 3. for each component, sum the soft assignment weights multiplied by the sample direction
// 4. get the length of (3) divide by the sum over samples of 2

using Float = float;
using Vector3 = Vec3f;

YBI_GPU bool SoftAssignment(const VMMDistribution &distribution,
                            const Vector3 &sampleDirection,
                            Float *softAssignment,
                            Float &total)
{
    Float V = 0.f;

    for (uint32_t componentIndex = 0; componentIndex < distribution.numComponents;
         componentIndex++)
    {
        Float v = distribution.weights[componentIndex] *
                  EvalVMF(sampleDirection,
                          distribution.meanDirection[componentIndex],
                          distribution.kappa[componentIndex]);

        V += v;

        softAssignment[componentIndex] = v;
    }

    const float divisor = V == 0.f ? 0.f : 1.f / V;
    total = V;
    for (uint32_t k = 0; k < distribution.numComponents; ++k)
    {
        softAssignment[k] *= divisor;
    }
    return V > 0.f;
}

struct GuidingParams
{
    const Vec3f *sampleDirection;
    const float *sampleRadianceEstimate;
    const float *samplePdf;
};

struct SampleData
{
    Vec3f sampleDirection;
    float radianceEstimate;
    float pdf;
};

struct SufficientStatistics
{
    // No default member initializers: the CUDA backend allocates this in shared memory.
    Vec3f sumWeightedDirections[VMM_MAX_NUM_COMPONENTS];
    float sumWeightedStats[VMM_MAX_NUM_COMPONENTS];

    // May be fractional after decaying previous observations.
    float numSamples;
};

struct VMMUpdateWorkspace
{
    SufficientStatistics currentStats;
    float logLikelihood;
    float sumWeights;
    float normalization;
};

struct VMMMapConfiguration
{
    float weightPrior{0.01f};
    float meanCosinePrior{0.f};
    float meanCosinePriorStrength{0.2f};
    float maxMeanCosine{1.f - 1.f / 32000.f};
};

struct VMMFitConfiguration
{
    VMMMapConfiguration map;
    uint32_t maxEMIterations{100};
    float convergenceThreshold{0.0025f};
    bool printLogLikelihood{false};
};

struct VMMUpdateRange
{
    uint32_t sampleOffset;
    uint32_t numSamples;
    // Lifetime count including this batch, before any decay of previous statistics.
    float overallNumSamples;
};

YBI_GPU SampleData LoadSampleData(const GuidingParams &params, uint32_t sampleIndex)
{
    SampleData result;
    result.sampleDirection = params.sampleDirection[sampleIndex];
    result.radianceEstimate = params.sampleRadianceEstimate[sampleIndex];
    result.pdf = params.samplePdf[sampleIndex];

    return result;
}

// Statistics must be normalized so that the assignment masses across components
// sum to numSamples for each batch, as in OpenPGL's SufficientStatistics::normalize.
// The caller schedules components and completes reductions before calling this.
// The combined current and previous sample count must be positive.
// overallNumSamples includes the current batch and previous samples before decay.
YBI_CPU_GPU void VMMMapStateUpdate(VMMDistribution &distribution,
                                   uint32_t componentIndex,
                                   const SufficientStatistics &currentStats,
                                   const SufficientStatistics &previousStats,
                                   float overallNumSamples,
                                   const VMMMapConfiguration &cfg = VMMMapConfiguration{})
{
    YBI_ASSERT(distribution.numComponents <= VMM_MAX_NUM_COMPONENTS);
    if (componentIndex >= distribution.numComponents)
    {
        return;
    }

    const float numSamples = currentStats.numSamples + previousStats.numSamples;
    YBI_ASSERT(numSamples > 0.f);

    const float newWeight =
        SafeDivide(cfg.weightPrior + currentStats.sumWeightedStats[componentIndex] +
                       previousStats.sumWeightedStats[componentIndex],
                   cfg.weightPrior * distribution.numComponents + numSamples);

    const Vec3f currentMean = SafeDivide(currentStats.sumWeightedDirections[componentIndex],
                                         currentStats.sumWeightedStats[componentIndex]);
    const Vec3f previousMean = SafeDivide(previousStats.sumWeightedDirections[componentIndex],
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

// All backend participants enter this function together. The caller provides
// shared workspace. The backend supplies worker indices, reductions, sync and reporting.
template <typename Backend>
YBI_GPU void UpdateVMM(const GuidingParams &params,
                       VMMDistribution &distribution,
                       SufficientStatistics &previousStats,
                       const VMMUpdateRange &range,
                       const VMMFitConfiguration &cfg,
                       VMMUpdateWorkspace &workspace,
                       Backend &backend)
{
    YBI_ASSERT(range.numSamples > 0);
    YBI_ASSERT(range.overallNumSamples >= float(range.numSamples));
    YBI_ASSERT(distribution.numComponents > 0 &&
               distribution.numComponents <= VMM_MAX_NUM_COMPONENTS);
    YBI_ASSERT(cfg.maxEMIterations > 0);
    YBI_ASSERT(cfg.convergenceThreshold >= 0.f);

    float previousLogLikelihood = 0.f;
    for (uint32_t iteration = 0; iteration < cfg.maxEMIterations; ++iteration)
    {
        // Finish reading the previous iteration before clearing shared results.
        backend.Sync();
        for (uint32_t k = backend.GetThreadIndex(); k < distribution.numComponents;
             k += backend.GetThreadStride())
        {
            workspace.currentStats.sumWeightedDirections[k] = Vec3f(0.f);
            workspace.currentStats.sumWeightedStats[k] = 0.f;
        }
        if (backend.GetThreadIndex() == 0)
        {
            workspace.currentStats.numSamples = float(range.numSamples);
            workspace.logLikelihood = 0.f;
            workspace.sumWeights = 0.f;
            workspace.normalization = 0.f;
        }
        backend.Sync();
        SufficientStatistics localStats{};
        float localLogLikelihood = 0.f;
        for (uint32_t sampleIndex = backend.GetThreadIndex(); sampleIndex < range.numSamples;
             sampleIndex += backend.GetThreadStride())
        {
            const SampleData sample = LoadSampleData(params, range.sampleOffset + sampleIndex);
            YBI_ASSERT(sample.pdf > 0.f);
            YBI_ASSERT(sample.radianceEstimate >= 0.f);
            const float sampleWeight = sample.radianceEstimate / sample.pdf;
            if (sampleWeight == 0.f)
            {
                continue;
            }
            float assignments[VMM_MAX_NUM_COMPONENTS];
            float pdf;
            if (!SoftAssignment(distribution, sample.sampleDirection, assignments, pdf))
            {
                continue;
            }
            localLogLikelihood += sampleWeight * logf(pdf);
            for (uint32_t k = 0; k < distribution.numComponents; ++k)
            {
                const float assignedWeight = assignments[k] * sampleWeight;
                localStats.sumWeightedDirections[k] += sample.sampleDirection * assignedWeight;
                localStats.sumWeightedStats[k] += assignedWeight;
            }
        }

        // Every worker reduces every component in the same order.
        for (uint32_t k = 0; k < distribution.numComponents; ++k)
        {
            backend.ReduceAdd(localStats.sumWeightedDirections[k],
                              workspace.currentStats.sumWeightedDirections[k]);
            backend.ReduceAdd(localStats.sumWeightedStats[k],
                              workspace.currentStats.sumWeightedStats[k]);
        }
        backend.ReduceAdd(localLogLikelihood, workspace.logLikelihood);
        backend.Sync();

        if (backend.GetThreadIndex() == 0)
        {
            float sum = 0.f;
            for (uint32_t k = 0; k < distribution.numComponents; ++k)
            {
                sum += workspace.currentStats.sumWeightedStats[k];
            }
            workspace.sumWeights = sum;
            workspace.normalization = SafeDivide(workspace.currentStats.numSamples, sum);
        }
        backend.Sync();
        // A batch with no positive assigned weight has no information to fit.
        if (workspace.sumWeights == 0.f)
        {
            return;
        }
        for (uint32_t k = backend.GetThreadIndex(); k < distribution.numComponents;
             k += backend.GetThreadStride())
        {
            workspace.currentStats.sumWeightedDirections[k] *= workspace.normalization;
            workspace.currentStats.sumWeightedStats[k] *= workspace.normalization;
        }
        backend.Sync();
        for (uint32_t componentIndex = backend.GetThreadIndex();
             componentIndex < distribution.numComponents;
             componentIndex += backend.GetThreadStride())
        {
            VMMMapStateUpdate(distribution,
                              componentIndex,
                              workspace.currentStats,
                              previousStats,
                              range.overallNumSamples,
                              cfg.map);
        }
        backend.Sync();

        const float logLikelihood = workspace.logLikelihood;
        backend.ReportIteration(
            iteration, logLikelihood, previousLogLikelihood, cfg.printLogLikelihood);
        if (iteration > 0 && Abs(logLikelihood - previousLogLikelihood) <=
                                 cfg.convergenceThreshold * Abs(previousLogLikelihood))
        {
            break;
        }
        previousLogLikelihood = logLikelihood;
    }
    // Commit only the final E-step statistics, not every EM iteration.
    for (uint32_t k = backend.GetThreadIndex(); k < distribution.numComponents;
         k += backend.GetThreadStride())
    {
        previousStats.sumWeightedDirections[k] += workspace.currentStats.sumWeightedDirections[k];
        previousStats.sumWeightedStats[k] += workspace.currentStats.sumWeightedStats[k];
    }
    if (backend.GetThreadIndex() == 0)
    {
        previousStats.numSamples += workspace.currentStats.numSamples;
    }
    backend.Sync();
}

} // namespace ybi
