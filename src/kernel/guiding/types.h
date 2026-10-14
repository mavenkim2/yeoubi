#pragma once

#include "core/lane.h"

namespace ybi
{

#define VMM_MAX_NUM_COMPONENTS 16
#define VMM_BLOCK_SIZE 256

constexpr uint32_t VMFSplitTableSize = 256;

struct ExclusiveScanResult
{
    uint32_t prefix;
    uint32_t total;
};

struct GuidingParams
{
    Vec3f *sampleDirection;
    float *sampleIncidentRadianceEstimate;
    float *samplePdf;
    uint64_t vmfSplitTexture;
    float vmfSplitMaxInvSqrtK;
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
    Vector3<Float> direction;
    Vector3<Float> position;
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

// Kd tree
enum class PathGuidingKdTreeStage
{
    PositionStatistics,
    Partition,
    Signature,
    Promotion,
};

struct PathGuidingKdTreeNode
{
    float splitLocation;
    int splitAxis;

    Vec3f boundsMin;
    Vec3f boundsMax;

    int64_t meanX;
    int64_t meanY;
    int64_t meanZ;

    int64_t varianceX;
    int64_t varianceY;
    int64_t varianceZ;
    uint32_t numSamples;

    int guidingCellIndex;
    int guidingNodeIndex = -1;

    YBI_CPU_GPU bool NeedsSplit() const
    {
        return splitAxis == -1;
    }

    YBI_CPU_GPU bool IsGuidingCell() const
    {
        return guidingCellIndex >= 0;
    }

    YBI_CPU_GPU bool IsLookahead() const
    {
        return !IsGuidingCell() && guidingNodeIndex >= 0;
    }

    YBI_CPU_GPU bool MatchesStage(PathGuidingKdTreeStage stage) const
    {
        switch (stage)
        {
        case PathGuidingKdTreeStage::PositionStatistics: return NeedsSplit();
        case PathGuidingKdTreeStage::Partition: return !NeedsSplit();
        case PathGuidingKdTreeStage::Signature: return IsGuidingCell() || IsLookahead();
        case PathGuidingKdTreeStage::Promotion: return IsLookahead();
        }
        return false;
    }
};

struct PathGuidingKdTreeWorkItem
{
    uint32_t nodeIndex;
    uint32_t sampleOffset;
    uint32_t sampleCount;

    // NOTE: stores count after PartitionCount, then the exclusive scan offset
    uint32_t leftOffset;
};

struct PathGuidingKdTreeNodeItem
{
    uint32_t nodeIndex;

    uint32_t sampleOffset;
    uint32_t sampleCount;
    uint32_t numLeft;

    uint32_t workItemOffset;
    uint32_t workItemCount;

    uint32_t treeDepth;
    uint32_t lookaheadDepth;
};

struct PathGuidingKdTreeParams
{
    PathGuidingKdTreeNode *nodes;

    PathGuidingKdTreeNodeItem *nodeItems;
    PathGuidingKdTreeWorkItem *workItems;

    uint32_t numNodeItems;
    uint32_t numWorkItems;

    SampleData *sampleData;
    uint32_t *sampleIndices;
};

} // namespace ybi
