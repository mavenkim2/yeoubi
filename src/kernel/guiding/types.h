#pragma once

#include "core/lane.h"

namespace ybi
{

#define VMM_MAX_NUM_COMPONENTS 16

struct GuidingParams
{
    Vec3f *sampleDirection;
    float *sampleIncidentRadianceEstimate;
    float *samplePdf;
};

struct VMMDistribution
{
    Vec3f meanDirection[VMM_MAX_NUM_COMPONENTS];
    // NOTE: must sum to 1
    float weights[VMM_MAX_NUM_COMPONENTS];
    float kappa[VMM_MAX_NUM_COMPONENTS];

    uint32_t numComponents;
};

template <typename Float>
struct VMMSampleData
{
    Vector3<Float> sampleDirection;
    Float incidentRadianceEstimate;
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
    float partialWeightScale;
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
    float splittingThreshold{0.75f};
    float mergingThreshold{0.00625f};
    bool printLogLikelihood{false};
};

struct VMMUpdateRange
{
    uint32_t sampleOffset;
    uint32_t numSamples;
    // Lifetime count including this batch, before any decay of previous statistics.
    float overallNumSamples;
};

template <typename Float>
struct VMMSplitStatistics
{
    Float chiSquareEstimates[VMM_MAX_NUM_COMPONENTS];
    Float sumWeights[VMM_MAX_NUM_COMPONENTS];
    Vector3<Float> sumWeightedCovariances[VMM_MAX_NUM_COMPONENTS];
};

struct VMMSplitWorkspace
{
    VMMSplitStatistics<float> statistics;
    Vec3f tangents[VMM_MAX_NUM_COMPONENTS];
    Vec3f bitangents[VMM_MAX_NUM_COMPONENTS];
    uint32_t numSplits;
    uint32_t componentMask;
};

struct VMMMergeWorkspace
{
    float costs[VMM_MAX_NUM_COMPONENTS * (VMM_MAX_NUM_COMPONENTS - 1) / 2];
    uint32_t numMerges;
    uint32_t componentMask;
    bool merged;
};

struct VMMMergeCandidate
{
    Vec3f meanDirection;
    float weight;
    float kappa;
};

} // namespace ybi
