#pragma once

#include "device/device.h"
#include "kernel/guiding/types.h"

namespace ybi
{

class PathGuidingSystem
{
    Device *device;
    PathGuidingSplitDecision splitDecision;
    const uint32_t maxLookaheadDepth = 6;
    DeviceBuffer<PathGuidingKdTreeParams> kernelParams;
    DeviceBuffer<uint32_t> activeNodeItems;
    DeviceBuffer<PathGuidingKdTreeNode> nodes;
    DeviceBuffer<PathGuidingKdTreeNodeItem> nodeItems;
    DeviceBuffer<uint32_t> nodeItemIndices;
    DeviceBuffer<uint32_t> nextNodeItemIndices;
    DeviceBuffer<uint32_t> workItemLeftOffsets;
    DeviceBuffer<PathGuidingKdTreePositionStatistics> positionStatistics;
    DeviceBuffer<PathGuidingKdTreePositionStatistics> newPositionStatistics;
    DeviceBuffer<PathGuidingKdTreeLookaheadStatistics> lookaheadStatistics;
    DeviceBuffer<PathGuidingKdTreeLookaheadStatistics> newLookaheadStatistics;
    DeviceBuffer<uint32_t> guidingCellSplitMask;
    DeviceBuffer<VMMDistribution> vmmDistributions;
    DeviceBuffer<SufficientStatistics> vmmStatistics;
    uint32_t numNodes = 1;
    uint32_t numVMMs = 0;
    uint32_t generation = 0;

    std::unique_ptr<DeviceKernel> generateNodes;
    std::unique_ptr<DeviceKernel> generateNodeItems;
    std::unique_ptr<DeviceKernel> generateWorkItems;
    std::unique_ptr<DeviceKernel> calculateNodeBounds;
    std::unique_ptr<DeviceKernel> calculateNodePositionStatistics;
    std::unique_ptr<DeviceKernel> generateReductionWorkItems;
    std::unique_ptr<DeviceKernel> reduceNodePositionStatistics;
    std::unique_ptr<DeviceKernel> accumulateNodePositionStatistics;
    std::unique_ptr<DeviceKernel> selectNodeSplits;
    std::unique_ptr<DeviceKernel> partitionCount;
    std::unique_ptr<DeviceKernel> partitionExclusiveScan;
    std::unique_ptr<DeviceKernel> partitionSamples;
    std::unique_ptr<DeviceKernel> updateLookaheadCell;
    std::unique_ptr<DeviceKernel> reduceLookaheadStatistics;
    std::unique_ptr<DeviceKernel> accumulateLookaheadCells;
    std::unique_ptr<DeviceKernel> shouldSplitGuidingCell;
    std::unique_ptr<DeviceKernel> promoteGuidingCells;
    std::unique_ptr<DeviceKernel> applyParallaxShift;

    template <typename T>
    void ReserveBuffer(DeviceBuffer<T> &buffer, size_t count, size_t preserveCount = 0);
    void UploadParams(PathGuidingKdTreeParams &params);
    void ReadBackParams(PathGuidingKdTreeParams &params);
    void Launch(DeviceKernel &kernel, uint32_t gridX, uint32_t blockX);
    void GenerateWorkItems(PathGuidingKdTreeParams &params, PathGuidingKdTreeStage stage);
    void ReserveFrontierBuffers(PathGuidingKdTreeParams &params);
    void ReservePositionStatisticsBuffers(PathGuidingKdTreeParams &params);
    void ReserveNodeItemBuffers(PathGuidingKdTreeParams &params);
    void ReserveLookaheadStatisticsBuffers(PathGuidingKdTreeParams &params);
    void ReserveVMMDistributions(PathGuidingKdTreeParams &params);
    template <typename T>
    void ReduceNodeStatistics(PathGuidingKdTreeParams &params,
                              DeviceKernel &kernel,
                              T *&statistics,
                              T *&newStatistics);
    void CalculatePositionStatistics(PathGuidingKdTreeParams &params);
    void SelectNodeSplits(const PathGuidingKdTreeParams &params);
    void PartitionSamples(PathGuidingKdTreeParams &params);
    void CreateNodeItems(PathGuidingKdTreeParams &params);
    void ResetGuidingCellSplits(PathGuidingKdTreeParams &params);
    void UpdateLookaheadCells(PathGuidingKdTreeParams &params);
    void SelectGuidingCellSplits(PathGuidingKdTreeParams &params);
    void PromoteGuidingCells(PathGuidingKdTreeParams &params);
    void UpdateSampleCountSplits(PathGuidingKdTreeParams &params);
    void UpdateLookaheadSplits(PathGuidingKdTreeParams &params);
    void BeginUpdate(PathGuidingKdTreeParams &params);
    void AdvanceFrontier(PathGuidingKdTreeParams &params);
    void FinishUpdate(PathGuidingKdTreeParams &params);

public:
    PathGuidingSystem(
        Device *device,
        DeviceKind deviceKind,
        PathGuidingSplitDecision splitDecision = PathGuidingSplitDecision::SampleCount);
    ~PathGuidingSystem();
    PathGuidingSystem(const PathGuidingSystem &) = delete;
    PathGuidingSystem &operator=(const PathGuidingSystem &) = delete;

    void UpdateKDTree(PathGuidingKdTreeParams &params);
    void UpdateVMMs(PathGuidingKdTreeParams &params);
};

} // namespace ybi
