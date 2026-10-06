#pragma once

#include "core/vec3f.h"
#include "device/cuda_compat.h"
#include "kernel/guiding/vmf.h"

#if !defined(__CUDACC__) && \
    (defined(__SSE__) || defined(_M_X64) || (defined(_M_IX86_FP) && _M_IX86_FP >= 1))
#include "core/simd_float4.h"
#endif

namespace ybi
{

// 1. calculate fluence estimate, which is just 1/N * sum of mcEstimate
// 2. calculate softassignment weight per sample per vmm component. multiply by sample weight
// 3. for each component, sum the soft assignment weights multiplied by the sample direction
// 4. get the length of (3) divide by the sum over samples of 2

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

struct GuidingParams
{
    const Vec3f *sampleDirection;
    const float *sampleRadianceEstimate;
    const float *samplePdf;
};

template <typename Float>
struct VMMSampleData
{
    Vector3<Float> sampleDirection;
    Float radianceEstimate;
    Float pdf;
};

using SampleData = VMMSampleData<float>;

template <typename Float>
struct VMMStatistics
{
    // No default member initializers: the CUDA backend allocates this in shared memory.
    Vector3<Float> sumWeightedDirections[VMM_MAX_NUM_COMPONENTS];
    Float sumWeightedStats[VMM_MAX_NUM_COMPONENTS];

    // May be fractional after decaying previous observations.
    float numSamples;
};

using SufficientStatistics = VMMStatistics<float>;

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

template <typename Sample>
YBI_GPU Sample LoadSampleData(const GuidingParams &params,
                              uint32_t sampleIndex,
                              uint32_t numValidSamples);

template <>
YBI_GPU SampleData LoadSampleData<SampleData>(const GuidingParams &params,
                                              uint32_t sampleIndex,
                                              uint32_t numValidSamples)
{
    YBI_ASSERT(numValidSamples == 1);
    return {params.sampleDirection[sampleIndex],
            params.sampleRadianceEstimate[sampleIndex],
            params.samplePdf[sampleIndex]};
}

#if !defined(__CUDACC__) && \
    (defined(__SSE__) || defined(_M_X64) || (defined(_M_IX86_FP) && _M_IX86_FP >= 1))
template <>
inline VMMSampleData<Float4> LoadSampleData<VMMSampleData<Float4>>(const GuidingParams &params,
                                                                   uint32_t sampleIndex,
                                                                   uint32_t numValidSamples)
{
    YBI_ASSERT(numValidSamples > 0 && numValidSamples <= Float4::Width);
    // Each range needs storage through sampleOffset + RoundUp(numSamples, 4).
    // Load full packets from padded arrays, including the final partial packet.
    VMMSampleData<Float4> sample{LoadUnaligned(params.sampleDirection + sampleIndex),
                                 LoadUnaligned(params.sampleRadianceEstimate + sampleIndex),
                                 LoadUnaligned(params.samplePdf + sampleIndex)};
    if (numValidSamples < Float4::Width)
    {
        const Mask4 valid{
            _mm_cmplt_ps(_mm_setr_ps(0.f, 1.f, 2.f, 3.f), _mm_set1_ps(float(numValidSamples)))};
        sample.sampleDirection.x = Select(valid, sample.sampleDirection.x, Float4(0.f));
        sample.sampleDirection.y = Select(valid, sample.sampleDirection.y, Float4(0.f));
        sample.sampleDirection.z = Select(valid, sample.sampleDirection.z, Float4(0.f));
        sample.radianceEstimate = Select(valid, sample.radianceEstimate, Float4(0.f));
        sample.pdf = Select(valid, sample.pdf, Float4(1.f));
    }
    return sample;
}
#endif

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

template <typename Backend>
YBI_GPU void BeginVMMIteration(uint32_t numSamples,
                               uint32_t numComponents,
                               VMMUpdateWorkspace &workspace,
                               Backend &backend)
{
    backend.Sync();
    for (uint32_t k = backend.GetThreadIndex(); k < numComponents; k += backend.GetThreadStride())
    {
        workspace.currentStats.sumWeightedDirections[k] = Vec3f(0.f);
        workspace.currentStats.sumWeightedStats[k] = 0.f;
    }
    if (backend.GetThreadIndex() == 0)
    {
        workspace.currentStats.numSamples = float(numSamples);
        workspace.logLikelihood = 0.f;
        workspace.sumWeights = 0.f;
        workspace.normalization = 0.f;
    }
    backend.Sync();
}

template <typename Backend>
YBI_GPU void ReduceVMMStatistics(const typename Backend::VMMStatistics &localStats,
                                 typename Backend::Float localLogLikelihood,
                                 uint32_t numComponents,
                                 VMMUpdateWorkspace &workspace,
                                 Backend &backend)
{
    for (uint32_t k = 0; k < numComponents; ++k)
    {
        backend.ReduceAdd(localStats.sumWeightedDirections[k],
                          workspace.currentStats.sumWeightedDirections[k]);
        backend.ReduceAdd(localStats.sumWeightedStats[k],
                          workspace.currentStats.sumWeightedStats[k]);
    }
    backend.ReduceAdd(localLogLikelihood, workspace.logLikelihood);
    backend.Sync();
}

template <typename Backend>
YBI_GPU bool
NormalizeVMMStatistics(uint32_t numComponents, VMMUpdateWorkspace &workspace, Backend &backend)
{
    if (backend.GetThreadIndex() == 0)
    {
        float sum = 0.f;
        for (uint32_t k = 0; k < numComponents; ++k)
        {
            sum += workspace.currentStats.sumWeightedStats[k];
        }
        workspace.sumWeights = sum;
        workspace.normalization = SafeDivide(workspace.currentStats.numSamples, sum);
    }
    backend.Sync();
    if (workspace.sumWeights == 0.f)
    {
        return false;
    }
    for (uint32_t k = backend.GetThreadIndex(); k < numComponents; k += backend.GetThreadStride())
    {
        workspace.currentStats.sumWeightedDirections[k] *= workspace.normalization;
        workspace.currentStats.sumWeightedStats[k] *= workspace.normalization;
    }
    backend.Sync();
    return true;
}

template <typename Backend>
YBI_GPU void UpdateVMMComponents(VMMDistribution &distribution,
                                 const SufficientStatistics &currentStats,
                                 const SufficientStatistics &previousStats,
                                 float overallNumSamples,
                                 const VMMMapConfiguration &cfg,
                                 Backend &backend)
{
    for (uint32_t componentIndex = backend.GetThreadIndex();
         componentIndex < distribution.numComponents;
         componentIndex += backend.GetThreadStride())
    {
        VMMMapStateUpdate(
            distribution, componentIndex, currentStats, previousStats, overallNumSamples, cfg);
    }
    backend.Sync();
}

template <typename Backend>
YBI_GPU void CommitVMMStatistics(SufficientStatistics &previousStats,
                                 const SufficientStatistics &currentStats,
                                 uint32_t numComponents,
                                 Backend &backend)
{
    for (uint32_t k = backend.GetThreadIndex(); k < numComponents; k += backend.GetThreadStride())
    {
        previousStats.sumWeightedDirections[k] += currentStats.sumWeightedDirections[k];
        previousStats.sumWeightedStats[k] += currentStats.sumWeightedStats[k];
    }
    if (backend.GetThreadIndex() == 0)
    {
        previousStats.numSamples += currentStats.numSamples;
    }
    backend.Sync();
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
        BeginVMMIteration(range.numSamples, distribution.numComponents, workspace, backend);
        VMMStatistics localStats{};
        Float localLogLikelihood(0.f);
        for (uint32_t sampleIndex = backend.GetThreadIndex() * laneWidth;
             sampleIndex < range.numSamples;
             sampleIndex += backend.GetThreadStride() * laneWidth)
        {
            const SampleData sample =
                LoadSampleData<SampleData>(params,
                                           range.sampleOffset + sampleIndex,
                                           Min(laneWidth, range.numSamples - sampleIndex));
            YBI_ASSERT(All(sample.pdf > Float(0.f)));
            YBI_ASSERT(All(sample.radianceEstimate >= Float(0.f)));
            const Float sampleWeight = sample.radianceEstimate / sample.pdf;
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

        ReduceVMMStatistics(
            localStats, localLogLikelihood, distribution.numComponents, workspace, backend);
        // A batch with no positive assigned weight has no information to fit.
        if (!NormalizeVMMStatistics(distribution.numComponents, workspace, backend))
        {
            return;
        }
        UpdateVMMComponents(distribution,
                            workspace.currentStats,
                            previousStats,
                            range.overallNumSamples,
                            cfg.map,
                            backend);

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
    CommitVMMStatistics(
        previousStats, workspace.currentStats, distribution.numComponents, backend);
}

} // namespace ybi
