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

struct VMMDistribution
{
    Vec3f meanDirection[VMM_MAX_NUM_COMPONENTS];
    // NOTE: must sum to 1
    float weights[VMM_MAX_NUM_COMPONENTS];
    float kappa[VMM_MAX_NUM_COMPONENTS];
    float distance[VMM_MAX_NUM_COMPONENTS];

    Vec3f pivotPosition;
    uint32_t numComponents;
    bool initialized;
    bool split;
};

struct PathGuidingTrainingSample
{
    Vec3f position;
    uint32_t octahedralDirection;
    float incidentRadiance;
    float pdf;
    float distance;
};

template <typename Float>
struct VMMSampleData
{
    Vector3<Float> direction;
    Vector3<Float> position;
    Float incidentRadianceEstimate;
    Float pdf;
    Float distance;
};

using SampleData = VMMSampleData<float>;

template <typename Float>
struct VMMStatistics
{
    // No default member initializers: the CUDA backend allocates this in shared memory.
    Vector3<Float> sumWeightedDirections[VMM_MAX_NUM_COMPONENTS];
    Float sumWeightedStats[VMM_MAX_NUM_COMPONENTS];
    Float sumDistanceWeights[VMM_MAX_NUM_COMPONENTS];

    // May be fractional after decaying previous observations.
    float numSamples;
    float sumWeights;
    float overallNumSamples;
    uint64_t numSamplesAfterLastSplit;
    uint64_t numSamplesAfterLastMerge;
};

using SufficientStatistics = VMMStatistics<float>;

struct VMMUpdateWorkspace
{
    SufficientStatistics currentStats;
    float logLikelihood;
    float normalization;
    float partialWeightScale;
};

struct VMMComponentDistanceWorkspace
{
    float batchSumWeights[VMM_MAX_NUM_COMPONENTS];
    float batchSumWeightedInverseDistances[VMM_MAX_NUM_COMPONENTS];
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
    uint32_t initialNumComponents{VMM_MAX_NUM_COMPONENTS / 2};
    float initialKappa{0.5f};
    bool useSplitAndMerge{true};
    bool partialReFit{true};
    uint32_t minSamplesForSplitting{4000};
    uint32_t minSamplesForPartialRefitting{4000};
    uint32_t minSamplesForMerging{8000};
    float decayOnSpatialSplit{0.25f};
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
    Float numSamples[VMM_MAX_NUM_COMPONENTS];
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
enum class PathGuidingSplitDecision
{
    SampleCount,
    Lookahead,
};

enum class PathGuidingKdTreeStage
{
    PositionStatistics,
    Partition,
    Lookahead,
    Promotion,
};

struct LookaheadCell
{
    float radianceFirstMoment;
    float radianceSecondMoment;
    uint32_t numSamples;
    Vec3f radianceWeightedDirectionSum;
};

struct PathGuidingKdTreePositionStatistics
{
    Vec3f mean;
    Vec3f variance;
    uint32_t numSamples;

    YBI_CPU_GPU PathGuidingKdTreePositionStatistics
    operator+(const PathGuidingKdTreePositionStatistics &other) const
    {
        if (numSamples == 0)
            return other;
        if (other.numSamples == 0)
            return *this;
        PathGuidingKdTreePositionStatistics result = {};
        result.numSamples = numSamples + other.numSamples;
        const float weightB = float(other.numSamples) / float(result.numSamples);
        const float weightA = float(numSamples) / float(result.numSamples);
        const Vec3f delta = other.mean - mean;
        result.mean = mean + delta * weightB;
        result.variance =
            variance * weightA + other.variance * weightB + delta * delta * (weightA * weightB);
        return result;
    }
};

struct PathGuidingKdTreeLookaheadStatistics
{
    float meanRadiance;
    float meanSquaredRadiance;
    Vec3f meanRadianceWeightedDirection;
    uint32_t numSamples;

    YBI_CPU_GPU PathGuidingKdTreeLookaheadStatistics
    operator+(const PathGuidingKdTreeLookaheadStatistics &other) const
    {
        if (numSamples == 0)
            return other;
        if (other.numSamples == 0)
            return *this;
        PathGuidingKdTreeLookaheadStatistics result = {};
        result.numSamples = numSamples + other.numSamples;
        const float weight = float(other.numSamples) / float(result.numSamples);
        result.meanRadiance = meanRadiance + (other.meanRadiance - meanRadiance) * weight;
        result.meanSquaredRadiance =
            meanSquaredRadiance + (other.meanSquaredRadiance - meanSquaredRadiance) * weight;
        result.meanRadianceWeightedDirection =
            meanRadianceWeightedDirection +
            (other.meanRadianceWeightedDirection - meanRadianceWeightedDirection) * weight;
        return result;
    }
};

struct PathGuidingKdTreeNode
{
    float splitLocation;
    int splitAxis = -1;
    int childOffset = -1;
    uint32_t treeDepth = 0;

    // Spatial region inherited from the parent, independent of sampled positions.
    Vec3f cellBoundsMin = Vec3f(-FLOAT_INFINITY);
    Vec3f cellBoundsMax = Vec3f(FLOAT_INFINITY);

    // Accumulated until a split plane is established; preserved on promotion.
    PathGuidingKdTreePositionStatistics positionStatistics = {};
    Vec3f positionBoundsMin = Vec3f(FLOAT_INFINITY);
    Vec3f positionBoundsMax = Vec3f(-FLOAT_INFINITY);
    uint32_t partitionGeneration = 0;

    int guidingCellIndex = -1;
    int ancestorGuidingCellIndex = -1;
    int parentIndex = -1;
    LookaheadCell lookaheadCell = {};

    YBI_CPU_GPU bool HasChildren() const
    {
        return childOffset >= 0;
    }

    YBI_CPU_GPU bool NeedsSplit() const
    {
        return splitAxis == -1;
    }

    YBI_CPU_GPU bool NeedsPartition(uint32_t generation) const
    {
        return !NeedsSplit() && partitionGeneration != generation;
    }

    YBI_CPU_GPU bool IsGuidingCell() const
    {
        return guidingCellIndex >= 0;
    }

    YBI_CPU_GPU int GetGuidingAncestorIndex(const PathGuidingKdTreeNode *nodes) const
    {
        int index = parentIndex;
        while (index >= 0 && !nodes[index].IsGuidingCell())
            index = nodes[index].parentIndex;
        return index;
    }

    YBI_CPU_GPU bool IsLookahead() const
    {
        return !IsGuidingCell() && ancestorGuidingCellIndex >= 0;
    }

    YBI_CPU_GPU bool MatchesStage(PathGuidingKdTreeStage stage, uint32_t generation) const
    {
        switch (stage)
        {
            case PathGuidingKdTreeStage::PositionStatistics:
                return NeedsSplit();
            case PathGuidingKdTreeStage::Partition:
                return NeedsPartition(generation);
            case PathGuidingKdTreeStage::Lookahead:
                return IsGuidingCell() || IsLookahead();
            case PathGuidingKdTreeStage::Promotion:
                return IsLookahead();
        }
        return false;
    }
};

struct PathGuidingKdTreeWorkItem
{
    uint32_t nodeIndex;
    uint32_t nodeItemIndex;
    uint32_t sampleOffset;
    uint32_t sampleCount;
};

struct PathGuidingKdTreeNodeItem
{
    uint32_t nodeIndex;

    uint32_t sampleOffset;
    uint32_t sampleCount;
    uint32_t numLeft;
    // Reset with the batch, preserved when this node is revisited after promotion.
    bool positionStatisticsAccumulated = false;

    uint32_t workItemOffset;
    uint32_t workItemCount;
    uint32_t reductionInputOffset;
    uint32_t reductionInputCount;
};

struct PathGuidingKdTreeParams
{
    PathGuidingSplitDecision splitDecision = PathGuidingSplitDecision::SampleCount;
    PathGuidingKdTreeStage stage = PathGuidingKdTreeStage::PositionStatistics;
    uint32_t maxLookaheadDepth = 0;
    PathGuidingKdTreeNode *nodes;
    uint32_t numNodes = 0;

    PathGuidingKdTreeNodeItem *nodeItems;
    uint32_t *nodeItemIndices;
    uint32_t numNodeItemIndices;
    uint32_t *nextNodeItemIndices;
    uint32_t numNextNodeItemIndices;
    uint32_t *workItemLeftOffsets;

    uint32_t numWorkItems;
    uint32_t maxWorkItemCount = 0;
    uint32_t *activeNodeItems;
    uint32_t numActiveNodeItems;
    uint32_t generation;
    PathGuidingKdTreePositionStatistics *positionStatistics;
    PathGuidingKdTreePositionStatistics *newPositionStatistics = nullptr;
    PathGuidingKdTreeLookaheadStatistics *lookaheadStatistics;
    PathGuidingKdTreeLookaheadStatistics *newLookaheadStatistics = nullptr;

    PathGuidingTrainingSample *sampleData;
    uint32_t *sampleIndices;
    uint32_t *newSampleIndices;
    uint32_t numSampleIndices;
    uint32_t *guidingCellSplitMask;
    uint32_t numGuidingCellSplits = 0;
    VMMDistribution *vmmDistributions;
    SufficientStatistics *vmmStatistics = nullptr;
    VMMSplitStatistics<float> *vmmSplitStatistics = nullptr;
    uint32_t numVMMs = 0;
    uint32_t *vmmUpdateNodeIndices = nullptr;
    uint32_t numVMMUpdateNodes = 0;
    uint32_t *vmmSplitComponentMasks = nullptr;
    VMMFitConfiguration vmmFitConfiguration;
    float meanRadianceSplitThreshold = 0.05f;

    uint64_t vmfSplitTexture;
    float vmfSplitMaxInvSqrtK;
};

} // namespace ybi
