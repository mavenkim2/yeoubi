#pragma once

#include "device/device.h"
#include "kernel/guiding/radiance_recorder.h"
#include "kernel/guiding/types.h"
#include "render/path_guiding/lut.h"

namespace ybi
{

class PathGuidingSystem
{
    Device *device;
    PathGuidingSplitDecision splitDecision;
    const uint32_t maxLookaheadDepth = 6;
    static constexpr uint32_t maxRadianceRecords =
        (256u * 1024u * 1024u) / (sizeof(RadianceRecord) + 2 * sizeof(int));
    DeviceBuffer<PathGuidingKdTreeParams> kernelParams;
    VMFSplitTable vmfSplitTable;
    RadianceRecorder radianceRecorder{};
    DeviceBuffer<RadianceRecorder> radianceRecorderParams;
    DeviceBuffer<RadianceRecord> radianceRecords;
    DeviceBuffer<int> radianceRecordFreeList;
    DeviceBuffer<int> radianceRecordRetiredList;
    DeviceBuffer<uint32_t> radianceRecordAllocationCursorBuffer;
    DeviceBuffer<uint32_t> radianceRecordRetiredCountBuffer;
    DeviceBuffer<PathGuidingTrainingSample> pathGuidingTrainingSamplesBuffer;
    DeviceBuffer<uint32_t> pathGuidingTrainingSampleCountBuffer;
    DeviceBuffer<PathGuidingTrainingSample> trainingSamples;
    DeviceBuffer<uint32_t> trainingSampleIndices;
    DeviceBuffer<uint32_t> trainingSampleScratchIndices;
    const size_t maxTrainingSampleBytes;
    size_t trainingSampleCount = 0;
    uint64_t droppedTrainingSampleCount = 0;
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
    DeviceBuffer<VMMSplitStatistics<float>> vmmSplitStatistics;
    DeviceBuffer<uint32_t> vmmUpdateNodeIndices;
    DeviceBuffer<uint32_t> vmmSplitComponentMasks;
    uint32_t numNodes = 1;
    uint32_t numVMMs = 0;
    uint32_t generation = 0;
    uint64_t renderedSamples = 0;
    uint64_t lastUpdateSample = 0;
    uint64_t nextUpdateSample = 1;
    uint32_t updateSampleCount = 1;
    uint32_t iteration = 0;

    struct StageTimer
    {
        const char *stage;
        std::unique_ptr<DeviceEventTimer> timer;
    };
    std::vector<StageTimer> stageTimers;
    bool timingUpdate = false;
    uint32_t updateFrontierCount = 0;

    std::unique_ptr<DeviceKernel> generateNodes;
    std::unique_ptr<DeviceKernel> initializeSampleIndices;
    std::unique_ptr<DeviceKernel> generateNodeItems;
    std::unique_ptr<DeviceKernel> copySampleCountDistributions;
    std::unique_ptr<DeviceKernel> generateWorkItems;
    std::unique_ptr<DeviceKernel> collectVMMUpdateNodes;
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
    std::unique_ptr<DeviceKernel> updateVMM;
    std::unique_ptr<DeviceKernel> partialUpdateVMM;
    std::unique_ptr<DeviceKernel> splitVMM;
    std::unique_ptr<DeviceKernel> mergeVMM;
    std::unique_ptr<DeviceKernel> updateComponentDistances;
    std::unique_ptr<DeviceKernel> applyParallaxShift;
    std::unique_ptr<DeviceKernel> reprojectSamples;
    std::unique_ptr<DeviceKernel> returnRetiredRecords;

    template <typename T>
    void ReserveBuffer(DeviceBuffer<T> &buffer, size_t count, size_t preserveCount = 0);
    void AllocateRadianceRecords(uint32_t count);
    void ReserveTrainingSampleBuffer();
    void DrainTrainingSamples();
    void PrintUpdateStatistics(size_t sampleCount) const;
    void UploadParams(PathGuidingKdTreeParams &params);
    void ReadBackParams(PathGuidingKdTreeParams &params);
    std::unique_ptr<DeviceEventTimer> StartStageTimer();
    void StopStageTimer(const char *stage, std::unique_ptr<DeviceEventTimer> timer);
    void PrintUpdateTimings(double wallMilliseconds, double lutMilliseconds);
    void Launch(DeviceKernel &kernel, uint32_t gridX, uint32_t blockX, const char *stage);
    void GenerateWorkItems(PathGuidingKdTreeParams &params, PathGuidingKdTreeStage stage);
    void ReserveFrontierBuffers(PathGuidingKdTreeParams &params);
    void ReservePositionStatisticsBuffers(PathGuidingKdTreeParams &params);
    void ReserveNodeItemBuffers(PathGuidingKdTreeParams &params);
    void ReserveLookaheadStatisticsBuffers(PathGuidingKdTreeParams &params);
    void ReserveVMMDistributions(PathGuidingKdTreeParams &params, uint32_t splitCount);
    void ReserveVMMTrainingBuffers(PathGuidingKdTreeParams &params);
    template <typename T>
    void ReduceNodeStatistics(PathGuidingKdTreeParams &params,
                              DeviceKernel &kernel,
                              T *&statistics,
                              T *&newStatistics,
                              const char *stage);
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
        PathGuidingSplitDecision splitDecision = PathGuidingSplitDecision::SampleCount,
        size_t maxTrainingSampleBytes = size_t(3) * 1024 * 1024 * 1024);
    ~PathGuidingSystem();
    PathGuidingSystem(const PathGuidingSystem &) = delete;
    PathGuidingSystem &operator=(const PathGuidingSystem &) = delete;

    void Update();
    void UpdateKDTree(PathGuidingKdTreeParams &params);
    void UpdateVMMs(PathGuidingKdTreeParams &params);

    RadianceRecorder BeginRecording(uint32_t recordCapacity);
    RadianceRecorder UpdateRadianceRecorder();
    void ClearTrainingSamples();

    DeviceBuffer<PathGuidingTrainingSample> GetTrainingSamples() const
    {
        return {trainingSamples.Data(), trainingSampleCount};
    }

    uint64_t GetDroppedTrainingSampleCount() const
    {
        return droppedTrainingSampleCount;
    }

    // Bind to KernelParams::pathGuiding on the same device after updating the field.
    const PathGuidingKdTreeParams *GetKernelParams() const
    {
        return kernelParams.Data();
    }
};

} // namespace ybi
