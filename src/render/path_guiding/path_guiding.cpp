#if defined(YBI_USE_PATH_GUIDING)

#if !defined(YBI_GUIDING_DISPATCH_KERNEL_PATH)
#error "YBI_USE_PATH_GUIDING requires YBI_GUIDING_DISPATCH_KERNEL_PATH"
#endif

#include "render/path_guiding/path_guiding.h"
#include "core/math.h"
#include "kernel/guiding/api.h"
#include "kernel/guiding/guiding.h"
#include "kernel/guiding/kd_tree.h"
#include "kernel/guiding/lookahead.h"
#include "kernel/guiding/radiance_recorder.h"
#include <chrono>
#include <cstdio>
#include <limits>
#include <map>
#include <string>
#include <utility>
#include <vector>

namespace ybi
{

PathGuidingSystem::PathGuidingSystem(Device *device,
                                     PathGuidingSplitDecision splitDecision,
                                     size_t maxTrainingSampleBytes)
    : device(device), splitDecision(splitDecision), maxTrainingSampleBytes(maxTrainingSampleBytes)
{
    YBI_ASSERT(device);
    DeviceKernelDesc desc;
    desc.filename = YBI_GUIDING_DISPATCH_KERNEL_PATH;
    desc.entryName = "ReturnRetiredRecordsKernel";
    returnRetiredRecords = device->CreateDeviceKernel(desc);
    YBI_ASSERT(returnRetiredRecords);
    radianceRecorderParams = device->Alloc<RadianceRecorder>(1);
    YBI_ASSERT(radianceRecorderParams);
    radianceRecordAllocationCursorBuffer = device->Alloc<uint32_t>(1);
    radianceRecordRetiredCountBuffer = device->Alloc<uint32_t>(1);
    pathGuidingTrainingSampleCountBuffer = device->Alloc<uint32_t>(1);
    radianceRecorder.allocationCursor = radianceRecordAllocationCursorBuffer.Data();
    radianceRecorder.retiredListCount = radianceRecordRetiredCountBuffer.Data();
    device->Memzero(radianceRecordAllocationCursorBuffer.Data(), sizeof(uint32_t));
    device->Memzero(radianceRecordRetiredCountBuffer.Data(), sizeof(uint32_t));
    device->Memzero(pathGuidingTrainingSampleCountBuffer.Data(), sizeof(uint32_t));
    desc.entryName = "InitializePathGuidingSampleIndicesKernel";
    initializeSampleIndices = device->CreateDeviceKernel(desc);
    YBI_ASSERT(initializeSampleIndices);
    desc.entryName = "CreateKdTreeWorkItemsKernel";
    generateWorkItems = device->CreateDeviceKernel(desc);
    desc.entryName = "CollectVMMUpdateNodesKernel";
    collectVMMUpdateNodes = device->CreateDeviceKernel(desc);
    desc.entryName = "CalculateNodeBoundsKernel";
    calculateNodeBounds = device->CreateDeviceKernel(desc);
    desc.entryName = "CalculateNodePositionStatisticsKernel";
    calculateNodePositionStatistics = device->CreateDeviceKernel(desc);
    desc.entryName = "CreateKdTreeReductionWorkItemsKernel";
    generateReductionWorkItems = device->CreateDeviceKernel(desc);
    desc.entryName = "ReduceNodePositionStatisticsKernel";
    reduceNodePositionStatistics = device->CreateDeviceKernel(desc);
    desc.entryName = "AccumulateNodePositionStatisticsKernel";
    accumulateNodePositionStatistics = device->CreateDeviceKernel(desc);
    desc.entryName = "SelectNodeSplitsKernel";
    selectNodeSplits = device->CreateDeviceKernel(desc);
    desc.entryName = "PartitionCountKernel";
    partitionCount = device->CreateDeviceKernel(desc);
    desc.entryName = "PartitionExclusiveScanKernel";
    partitionExclusiveScan = device->CreateDeviceKernel(desc);
    desc.entryName = "PartitionSamplesKernel";
    partitionSamples = device->CreateDeviceKernel(desc);
    desc.entryName = "CreateKdTreeNodesKernel";
    generateNodes = device->CreateDeviceKernel(desc);
    desc.entryName = "CreateKdTreeNodeItemsKernel";
    generateNodeItems = device->CreateDeviceKernel(desc);
    YBI_ASSERT(generateWorkItems && collectVMMUpdateNodes && calculateNodeBounds &&
               calculateNodePositionStatistics && generateReductionWorkItems &&
               reduceNodePositionStatistics && accumulateNodePositionStatistics &&
               selectNodeSplits && partitionCount && partitionExclusiveScan && partitionSamples &&
               generateNodes && generateNodeItems);

    desc.entryName = "UpdateVMMKernel";
    updateVMM = device->CreateDeviceKernel(desc);
    desc.entryName = "PartialUpdateVMMKernel";
    partialUpdateVMM = device->CreateDeviceKernel(desc);
    desc.entryName = "SplitVMMKernel";
    splitVMM = device->CreateDeviceKernel(desc);
    desc.entryName = "MergeVMMKernel";
    mergeVMM = device->CreateDeviceKernel(desc);
    desc.entryName = "UpdateComponentDistancesKernel";
    updateComponentDistances = device->CreateDeviceKernel(desc);
    desc.entryName = "ApplyParallaxShiftDispatchKernel";
    applyParallaxShift = device->CreateDeviceKernel(desc);
    desc.entryName = "ReprojectSamplesKernel";
    reprojectSamples = device->CreateDeviceKernel(desc);
    YBI_ASSERT(updateVMM && partialUpdateVMM && splitVMM && mergeVMM && updateComponentDistances &&
               applyParallaxShift && reprojectSamples);

    switch (splitDecision)
    {
        case PathGuidingSplitDecision::SampleCount:
            desc.entryName = "CopySampleCountDistributionsKernel";
            copySampleCountDistributions = device->CreateDeviceKernel(desc);
            YBI_ASSERT(copySampleCountDistributions);
            break;
        case PathGuidingSplitDecision::Lookahead:
            desc.entryName = "UpdateLookaheadCellKernel";
            updateLookaheadCell = device->CreateDeviceKernel(desc);
            desc.entryName = "ReduceLookaheadStatisticsKernel";
            reduceLookaheadStatistics = device->CreateDeviceKernel(desc);
            desc.entryName = "AccumulateLookaheadCellsKernel";
            accumulateLookaheadCells = device->CreateDeviceKernel(desc);
            desc.entryName = "ShouldSplitGuidingCellKernel";
            shouldSplitGuidingCell = device->CreateDeviceKernel(desc);
            desc.entryName = "PromoteGuidingCellsKernel";
            promoteGuidingCells = device->CreateDeviceKernel(desc);
            YBI_ASSERT(updateLookaheadCell && reduceLookaheadStatistics &&
                       accumulateLookaheadCells && shouldSplitGuidingCell && promoteGuidingCells);
            break;
    }
    kernelParams = device->Alloc<PathGuidingKdTreeParams>(1);
    YBI_ASSERT(kernelParams);
    PathGuidingKdTreeParams initialParams{};
    UploadParams(initialParams);

    PathGuidingKdTreeNode root{};
    root.guidingCellIndex = 0;
    VMMDistribution rootDistribution{};
    InitializeUniformVMM(rootDistribution);
    vmmDistributions = device->Alloc<VMMDistribution>(1);
    device->CopyToDevice(vmmDistributions.Data(), &rootDistribution, sizeof(rootDistribution));
    SufficientStatistics rootStatistics{};
    vmmStatistics = device->Alloc<SufficientStatistics>(1);
    device->CopyToDevice(vmmStatistics.Data(), &rootStatistics, sizeof(rootStatistics));
    VMMSplitStatistics<float> rootSplitStatistics{};
    vmmSplitStatistics = device->Alloc<VMMSplitStatistics<float>>(1);
    device->CopyToDevice(
        vmmSplitStatistics.Data(), &rootSplitStatistics, sizeof(rootSplitStatistics));
    numVMMs = 1;
    nodes = device->Alloc<PathGuidingKdTreeNode>(1);
    device->CopyToDevice(nodes.Data(), &root, sizeof(root));

    PathGuidingKdTreeNodeItem rootItem{};
    nodeItems = device->Alloc<PathGuidingKdTreeNodeItem>(1);
    device->CopyToDevice(nodeItems.Data(), &rootItem, sizeof(rootItem));
    uint32_t rootItemIndex = 0;
    nodeItemIndices = device->Alloc<uint32_t>(1);
    device->CopyToDevice(nodeItemIndices.Data(), &rootItemIndex, sizeof(rootItemIndex));
}

PathGuidingSystem::~PathGuidingSystem()
{
    ClearTrainingSamples();
    if (vmfSplitTable.texture.textureObject)
        device->DestroyTexture(vmfSplitTable.texture);
    if (trainingSampleIndices)
        device->Free(trainingSampleIndices);
    if (trainingSampleScratchIndices)
        device->Free(trainingSampleScratchIndices);
    if (radianceRecorderParams)
        device->Free(radianceRecorderParams);
    if (radianceRecords)
        device->Free(radianceRecords);
    if (radianceRecordFreeList)
        device->Free(radianceRecordFreeList);
    if (radianceRecordRetiredList)
        device->Free(radianceRecordRetiredList);
    if (radianceRecordAllocationCursorBuffer)
        device->Free(radianceRecordAllocationCursorBuffer);
    if (radianceRecordRetiredCountBuffer)
        device->Free(radianceRecordRetiredCountBuffer);
    if (pathGuidingTrainingSamplesBuffer)
        device->Free(pathGuidingTrainingSamplesBuffer);
    if (pathGuidingTrainingSampleCountBuffer)
        device->Free(pathGuidingTrainingSampleCountBuffer);
    if (kernelParams)
        device->Free(kernelParams);
    if (activeNodeItems)
        device->Free(activeNodeItems);
    if (nodes)
        device->Free(nodes);
    if (nodeItems)
        device->Free(nodeItems);
    if (nodeItemIndices)
        device->Free(nodeItemIndices);
    if (nextNodeItemIndices)
        device->Free(nextNodeItemIndices);
    if (workItemLeftOffsets)
        device->Free(workItemLeftOffsets);
    if (positionStatistics)
        device->Free(positionStatistics);
    if (newPositionStatistics)
        device->Free(newPositionStatistics);
    if (lookaheadStatistics)
        device->Free(lookaheadStatistics);
    if (newLookaheadStatistics)
        device->Free(newLookaheadStatistics);
    if (guidingCellSplitMask)
        device->Free(guidingCellSplitMask);
    if (vmmDistributions)
        device->Free(vmmDistributions);
    if (vmmStatistics)
        device->Free(vmmStatistics);
    if (vmmSplitStatistics)
        device->Free(vmmSplitStatistics);
    if (vmmUpdateNodeIndices)
        device->Free(vmmUpdateNodeIndices);
    if (vmmSplitComponentMasks)
        device->Free(vmmSplitComponentMasks);
}

template <typename T>
void PathGuidingSystem::ReserveBuffer(DeviceBuffer<T> &buffer, size_t count, size_t preserveCount)
{
    if (buffer.Count() >= count)
        return;
    DeviceBuffer<T> replacement = device->Alloc<T>(Max(count, buffer.Count() * 2));
    if (preserveCount != 0)
    {
        device->CopyDeviceToDeviceAsync(
            replacement.Data(), buffer.Data(), preserveCount * sizeof(T));
        device->SyncDeviceStream();
    }
    device->Free(buffer);
    buffer = replacement;
}

void PathGuidingSystem::UploadParams(PathGuidingKdTreeParams &params)
{
    device->CopyToDevice(kernelParams.Data(), &params, sizeof(params));
}

void PathGuidingSystem::ReadBackParams(PathGuidingKdTreeParams &params)
{
    device->SyncDeviceStream();
    device->CopyToHost(&params, kernelParams.Data(), sizeof(params));
}

std::unique_ptr<DeviceEventTimer> PathGuidingSystem::StartStageTimer()
{
    return timingUpdate ? device->StartEventTimer() : nullptr;
}

void PathGuidingSystem::StopStageTimer(const char *stage, std::unique_ptr<DeviceEventTimer> timer)
{
    if (timer)
    {
        device->StopEventTimer(*timer);
        stageTimers.push_back({stage, std::move(timer)});
    }
}

void PathGuidingSystem::Launch(DeviceKernel &kernel,
                                uint32_t gridX,
                                uint32_t blockX,
                                const char *stage)
{
    void *arguments[] = {kernelParams.KernelArg()};
    DeviceKernelLaunchDesc desc;
    desc.arguments = arguments;
    desc.gridX = gridX;
    desc.blockX = blockX;
    auto timer = StartStageTimer();
    device->ExecuteKernel(kernel, desc);
    StopStageTimer(stage, std::move(timer));
}

void PathGuidingSystem::AllocateRadianceRecords(uint32_t count)
{
    count = Min(count, maxRadianceRecords);
    const uint32_t oldCount = radianceRecorder.availableCount;
    if (count <= oldCount)
        return;

    DeviceBuffer<RadianceRecord> newRecords = device->Alloc<RadianceRecord>(count);
    DeviceBuffer<int> newFreeList = device->Alloc<int>(count);
    DeviceBuffer<int> newRetiredList = device->Alloc<int>(count);
    if (oldCount != 0)
    {
        device->CopyDeviceToDeviceAsync(
            newRecords.Data(), radianceRecords.Data(), size_t(oldCount) * sizeof(RadianceRecord));
        device->CopyDeviceToDeviceAsync(
            newFreeList.Data(), radianceRecordFreeList.Data(), size_t(oldCount) * sizeof(int));
        device->SyncDeviceStream();
    }

    std::vector<int> newIndices(count - oldCount);
    for (uint32_t i = oldCount; i < count; ++i)
        newIndices[i - oldCount] = int(i);
    device->CopyToDevice(
        newFreeList.Data() + oldCount, newIndices.data(), newIndices.size() * sizeof(int));

    device->Free(radianceRecords);
    device->Free(radianceRecordFreeList);
    device->Free(radianceRecordRetiredList);
    radianceRecords = newRecords;
    radianceRecordFreeList = newFreeList;
    radianceRecordRetiredList = newRetiredList;
    radianceRecorder.records = radianceRecords.Data();
    radianceRecorder.freeList = radianceRecordFreeList.Data();
    radianceRecorder.retiredList = radianceRecordRetiredList.Data();
    radianceRecorder.availableCount = count;
    radianceRecorder.retiredListMax = count;
}

void PathGuidingSystem::ReserveTrainingSampleBuffer()
{
    if (pathGuidingTrainingSamplesBuffer.Count() < radianceRecorder.availableCount)
    {
        device->SyncDeviceStream();
        if (pathGuidingTrainingSamplesBuffer)
            device->Free(pathGuidingTrainingSamplesBuffer);
        pathGuidingTrainingSamplesBuffer =
            device->Alloc<PathGuidingTrainingSample>(radianceRecorder.availableCount);
    }
    radianceRecorder.trainingSamples = pathGuidingTrainingSamplesBuffer.Data();
    radianceRecorder.trainingSampleCount = pathGuidingTrainingSampleCountBuffer.Data();
    radianceRecorder.trainingSamplesMax =
        static_cast<uint32_t>(pathGuidingTrainingSamplesBuffer.Count());
}

void PathGuidingSystem::DrainTrainingSamples()
{
    uint32_t writtenCount;
    device->CopyToHost(
        &writtenCount, pathGuidingTrainingSampleCountBuffer.Data(), sizeof(writtenCount));
    const uint32_t sampleCount =
        Min(writtenCount, uint32_t(pathGuidingTrainingSamplesBuffer.Count()));
    droppedTrainingSampleCount += writtenCount - sampleCount;

    const size_t maxSampleCount = maxTrainingSampleBytes / sizeof(PathGuidingTrainingSample);
    const size_t appendCount = Min<size_t>(sampleCount, maxSampleCount - trainingSampleCount);
    if (appendCount != 0)
    {
        const size_t requiredCount = trainingSampleCount + appendCount;
        if (trainingSamples.Count() < requiredCount)
        {
            const size_t capacity =
                Min(maxSampleCount, Max(trainingSampleCount * 2, requiredCount));
            DeviceBuffer<PathGuidingTrainingSample> replacement =
                device->Alloc<PathGuidingTrainingSample>(capacity);
            YBI_ASSERT(replacement);
            if (trainingSampleCount != 0)
            {
                device->CopyDeviceToDeviceAsync(replacement.Data(),
                                                trainingSamples.Data(),
                                                trainingSampleCount *
                                                    sizeof(PathGuidingTrainingSample));
                device->SyncDeviceStream();
            }
            device->Free(trainingSamples);
            trainingSamples = replacement;
        }
        device->CopyDeviceToDeviceAsync(trainingSamples.Data() + trainingSampleCount,
                                        pathGuidingTrainingSamplesBuffer.Data(),
                                        appendCount * sizeof(PathGuidingTrainingSample));
        trainingSampleCount = requiredCount;
    }
    droppedTrainingSampleCount += sampleCount - appendCount;
    device->Memzero(pathGuidingTrainingSampleCountBuffer.Data(), sizeof(uint32_t));
}

void PathGuidingSystem::ClearTrainingSamples()
{
    device->SyncDeviceStream();
    device->Free(trainingSamples);
    trainingSampleCount = 0;
    droppedTrainingSampleCount = 0;
    device->Memzero(pathGuidingTrainingSampleCountBuffer.Data(), sizeof(uint32_t));
}

RadianceRecorder PathGuidingSystem::BeginRecording(uint32_t recordCapacity)
{
    device->SyncDeviceStream();
    DrainTrainingSamples();
    AllocateRadianceRecords(recordCapacity);
    ReserveTrainingSampleBuffer();
    radianceRecorder.allocationCursor = radianceRecordAllocationCursorBuffer.Data();
    radianceRecorder.retiredListCount = radianceRecordRetiredCountBuffer.Data();
    device->Memzero(radianceRecordAllocationCursorBuffer.Data(), sizeof(uint32_t));
    device->Memzero(radianceRecordRetiredCountBuffer.Data(), sizeof(uint32_t));
    device->Memzero(pathGuidingTrainingSampleCountBuffer.Data(), sizeof(uint32_t));
    device->CopyToDevice(
        radianceRecorderParams.Data(), &radianceRecorder, sizeof(radianceRecorder));
    return radianceRecorder;
}

RadianceRecorder PathGuidingSystem::UpdateRadianceRecorder()
{
    YBI_ASSERT(radianceRecorder.allocationCursor && radianceRecorder.retiredListCount);
    YBI_ASSERT(radianceRecorder.availableCount == 0 || radianceRecorder.freeList);
    YBI_ASSERT(radianceRecorder.retiredListMax == 0 || radianceRecorder.retiredList);

    uint32_t allocationCount;
    device->SyncDeviceStream();
    device->CopyToHost(
        &allocationCount, radianceRecorder.allocationCursor, sizeof(allocationCount));
    DrainTrainingSamples();

    void *arguments[] = {radianceRecorderParams.KernelArg()};
    DeviceKernelLaunchDesc desc;
    desc.arguments = arguments;
    desc.gridX = 1;
    desc.blockX = VMM_BLOCK_SIZE;
    device->ExecuteKernel(*returnRetiredRecords, desc);

    if (allocationCount > radianceRecorder.availableCount &&
        radianceRecorder.availableCount < maxRadianceRecords)
    {
        device->SyncDeviceStream();
        const uint32_t capacity = uint32_t(Min<uint64_t>(
            maxRadianceRecords,
            Max<uint64_t>(allocationCount, uint64_t(radianceRecorder.availableCount) * 2)));
        AllocateRadianceRecords(capacity);
        ReserveTrainingSampleBuffer();
        device->CopyToDevice(
            radianceRecorderParams.Data(), &radianceRecorder, sizeof(radianceRecorder));
    }
    return radianceRecorder;
}

void PathGuidingSystem::GenerateWorkItems(PathGuidingKdTreeParams &params,
                                          PathGuidingKdTreeStage stage)
{
    params.stage = stage;
    UploadParams(params);
    Launch(*generateWorkItems, 1, VMM_BLOCK_SIZE, "kd.work_items");
    ReadBackParams(params);
}

void PathGuidingSystem::ReserveFrontierBuffers(PathGuidingKdTreeParams &params)
{
    ReserveBuffer(activeNodeItems, params.numNodeItemIndices);
    const size_t workItemCapacity =
        CeilDiv(params.numSampleIndices, VMM_KD_TREE_WORK_ITEM_SAMPLE_LIMIT) +
        size_t(params.numNodeItemIndices);
    ReserveBuffer(workItemLeftOffsets, workItemCapacity);
    params.activeNodeItems = activeNodeItems.Data();
    params.workItemLeftOffsets = workItemLeftOffsets.Data();
}

void PathGuidingSystem::ReservePositionStatisticsBuffers(PathGuidingKdTreeParams &params)
{
    ReserveBuffer(positionStatistics, params.numWorkItems);
    ReserveBuffer(newPositionStatistics, params.numWorkItems);
    params.positionStatistics = positionStatistics.Data();
    params.newPositionStatistics = newPositionStatistics.Data();
}

void PathGuidingSystem::ReserveNodeItemBuffers(PathGuidingKdTreeParams &params)
{
    const size_t nodeCapacity = size_t(params.numNodes) + size_t(2) * params.numNodeItemIndices;
    ReserveBuffer(nodes, nodeCapacity, params.numNodes);
    ReserveBuffer(nodeItems, nodeCapacity, params.numNodes);
    ReserveBuffer(nextNodeItemIndices,
                  size_t(2) * (size_t(params.numNodeItemIndices) + params.numGuidingCellSplits));
    params.nodes = nodes.Data();
    params.nodeItems = nodeItems.Data();
    params.nextNodeItemIndices = nextNodeItemIndices.Data();
}

void PathGuidingSystem::ReserveLookaheadStatisticsBuffers(PathGuidingKdTreeParams &params)
{
    ReserveBuffer(lookaheadStatistics, params.numWorkItems);
    ReserveBuffer(newLookaheadStatistics, params.numWorkItems);
    params.lookaheadStatistics = lookaheadStatistics.Data();
    params.newLookaheadStatistics = newLookaheadStatistics.Data();
}

void PathGuidingSystem::ReserveVMMDistributions(PathGuidingKdTreeParams &params,
                                                uint32_t splitCount)
{
    ReserveBuffer(vmmDistributions, size_t(params.numVMMs) + splitCount, params.numVMMs);
    ReserveBuffer(vmmStatistics, size_t(params.numVMMs) + splitCount, params.numVMMs);
    ReserveBuffer(vmmSplitStatistics, size_t(params.numVMMs) + splitCount, params.numVMMs);
    params.vmmDistributions = vmmDistributions.Data();
    params.vmmStatistics = vmmStatistics.Data();
    params.vmmSplitStatistics = vmmSplitStatistics.Data();
}

void PathGuidingSystem::ReserveVMMTrainingBuffers(PathGuidingKdTreeParams &params)
{
    ReserveBuffer(vmmSplitComponentMasks, params.numVMMs);
    params.vmmSplitComponentMasks = vmmSplitComponentMasks.Data();
}

template <typename T>
void PathGuidingSystem::ReduceNodeStatistics(PathGuidingKdTreeParams &params,
                                             DeviceKernel &kernel,
                                             T *&statistics,
                                             T *&newStatistics,
                                             const char *stage)
{
    uint32_t workItemBound = params.numWorkItems;
    for (uint32_t count = params.maxWorkItemCount; count > 1;
         count = CeilDiv(count, VMM_BLOCK_SIZE))
    {
        std::swap(statistics, newStatistics);
        UploadParams(params);
        Launch(*generateReductionWorkItems, 1, 1, "kd.reduction_setup");
        workItemBound = CeilDiv(workItemBound, VMM_BLOCK_SIZE) + params.numActiveNodeItems;
        Launch(kernel, workItemBound, VMM_BLOCK_SIZE, stage);
    }
}

void PathGuidingSystem::CalculatePositionStatistics(PathGuidingKdTreeParams &params)
{
    GenerateWorkItems(params, PathGuidingKdTreeStage::PositionStatistics);
    if (params.numWorkItems == 0)
        return;

    ReservePositionStatisticsBuffers(params);
    UploadParams(params);

    Launch(*calculateNodeBounds, params.numWorkItems, VMM_BLOCK_SIZE, "kd.bounds");
    Launch(*calculateNodePositionStatistics, params.numWorkItems, VMM_BLOCK_SIZE,
           "kd.position_stats");

    ReduceNodeStatistics(params,
                         *reduceNodePositionStatistics,
                         params.positionStatistics,
                         params.newPositionStatistics,
                         "kd.position_stats");
    Launch(*accumulateNodePositionStatistics,
           CeilDiv(params.numActiveNodeItems, VMM_BLOCK_SIZE),
           VMM_BLOCK_SIZE,
           "kd.position_stats");
}

void PathGuidingSystem::SelectNodeSplits(const PathGuidingKdTreeParams &params)
{
    if (params.numWorkItems == 0)
        return;

    Launch(*selectNodeSplits, CeilDiv(params.numActiveNodeItems, VMM_BLOCK_SIZE), VMM_BLOCK_SIZE,
           "kd.split_decision");
}

void PathGuidingSystem::PartitionSamples(PathGuidingKdTreeParams &params)
{
    GenerateWorkItems(params, PathGuidingKdTreeStage::Partition);
    if (params.numWorkItems == 0)
        return;

    YBI_ASSERT(params.newSampleIndices && params.newSampleIndices != params.sampleIndices);
    Launch(*partitionCount, params.numWorkItems, VMM_BLOCK_SIZE, "kd.partition");
    Launch(*partitionExclusiveScan, params.numActiveNodeItems, VMM_BLOCK_SIZE, "kd.partition");

    auto copyTimer = StartStageTimer();
    device->CopyDeviceToDeviceAsync(params.newSampleIndices,
                                    params.sampleIndices,
                                    size_t(params.numSampleIndices) * sizeof(uint32_t));
    StopStageTimer("kd.partition_copy", std::move(copyTimer));
    Launch(*partitionSamples, params.numWorkItems, VMM_BLOCK_SIZE, "kd.partition");
    std::swap(params.sampleIndices, params.newSampleIndices);
}

void PathGuidingSystem::CreateNodeItems(PathGuidingKdTreeParams &params)
{
    ReserveNodeItemBuffers(params);
    UploadParams(params);
    Launch(*generateNodes, 1, 1, "kd.nodes");
    Launch(*generateNodeItems, 1, 1, "kd.nodes");
    ReadBackParams(params);
}

void PathGuidingSystem::UpdateSampleCountSplits(PathGuidingKdTreeParams &params)
{
    CreateNodeItems(params);
    ReserveVMMDistributions(params, params.numNodeItemIndices);
    UploadParams(params);
    Launch(*copySampleCountDistributions, 1, 1, "kd.copy_distributions");
    ReadBackParams(params);
}

void PathGuidingSystem::ResetGuidingCellSplits(PathGuidingKdTreeParams &params)
{
    const size_t nodeCapacity = size_t(params.numNodes) + size_t(2) * params.numNodeItemIndices;
    const size_t maskWordCount = CeilDiv(nodeCapacity, size_t(32));
    ReserveBuffer(guidingCellSplitMask, maskWordCount);
    device->Memzero(guidingCellSplitMask.Data(), maskWordCount * sizeof(uint32_t));
    params.guidingCellSplitMask = guidingCellSplitMask.Data();
    params.numGuidingCellSplits = 0;
}

void PathGuidingSystem::UpdateLookaheadCells(PathGuidingKdTreeParams &params)
{
    GenerateWorkItems(params, PathGuidingKdTreeStage::Lookahead);
    if (params.numWorkItems == 0)
        return;

    ReserveLookaheadStatisticsBuffers(params);
    UploadParams(params);
    Launch(*updateLookaheadCell, params.numWorkItems, VMM_BLOCK_SIZE, "kd.lookahead_stats");

    ReduceNodeStatistics(params,
                         *reduceLookaheadStatistics,
                         params.lookaheadStatistics,
                         params.newLookaheadStatistics,
                         "kd.lookahead_stats");
    Launch(*accumulateLookaheadCells,
           CeilDiv(params.numActiveNodeItems, VMM_BLOCK_SIZE),
           VMM_BLOCK_SIZE,
           "kd.lookahead_stats");
}

void PathGuidingSystem::SelectGuidingCellSplits(PathGuidingKdTreeParams &params)
{
    if (params.numWorkItems == 0)
        return;

    Launch(*shouldSplitGuidingCell,
           CeilDiv(params.numActiveNodeItems, VMM_BLOCK_SIZE),
           VMM_BLOCK_SIZE,
           "kd.lookahead_decision");
    ReadBackParams(params);
}

void PathGuidingSystem::PromoteGuidingCells(PathGuidingKdTreeParams &params)
{
    if (params.numGuidingCellSplits == 0)
        return;

    ReserveVMMDistributions(params, params.numGuidingCellSplits);
    UploadParams(params);
    Launch(*promoteGuidingCells, 1, VMM_BLOCK_SIZE, "kd.promote");
    ReadBackParams(params);
}

void PathGuidingSystem::UpdateLookaheadSplits(PathGuidingKdTreeParams &params)
{
    ResetGuidingCellSplits(params);
    UpdateLookaheadCells(params);
    SelectGuidingCellSplits(params);
    CreateNodeItems(params);
    PromoteGuidingCells(params);
}

void PathGuidingSystem::BeginUpdate(PathGuidingKdTreeParams &params)
{
    params.generation = ++generation;
    params.splitDecision = splitDecision;
    params.maxLookaheadDepth = maxLookaheadDepth;
    params.numWorkItems = 0;
    params.numActiveNodeItems = 0;
    params.numNextNodeItemIndices = 0;
    params.numGuidingCellSplits = 0;
    YBI_ASSERT(params.numSampleIndices == 0 || (params.sampleData && params.sampleIndices));

    ReserveBuffer(nodeItems, numNodes);
    device->Memzero(nodeItems.Data(), size_t(numNodes) * sizeof(PathGuidingKdTreeNodeItem));
    PathGuidingKdTreeNodeItem rootItem{};
    rootItem.sampleCount = params.numSampleIndices;
    device->CopyToDevice(nodeItems.Data(), &rootItem, sizeof(rootItem));
    ReserveBuffer(nodeItemIndices, 1);
    uint32_t rootItemIndex = 0;
    device->CopyToDevice(nodeItemIndices.Data(), &rootItemIndex, sizeof(rootItemIndex));
    params.nodes = nodes.Data();
    params.numNodes = numNodes;
    params.nodeItems = nodeItems.Data();
    params.nodeItemIndices = nodeItemIndices.Data();
    params.numNodeItemIndices = 1;
    params.nextNodeItemIndices = nextNodeItemIndices.Data();
    params.guidingCellSplitMask = guidingCellSplitMask.Data();
    params.vmmDistributions = vmmDistributions.Data();
    params.vmmStatistics = vmmStatistics.Data();
    params.vmmSplitStatistics = vmmSplitStatistics.Data();
    params.numVMMs = numVMMs;
}

void PathGuidingSystem::AdvanceFrontier(PathGuidingKdTreeParams &params)
{
    params.numNodeItemIndices = params.numNextNodeItemIndices;
    params.numNextNodeItemIndices = 0;
    std::swap(nodeItemIndices, nextNodeItemIndices);
    params.nodeItemIndices = nodeItemIndices.Data();
    params.nextNodeItemIndices = nextNodeItemIndices.Data();
}

void PathGuidingSystem::FinishUpdate(PathGuidingKdTreeParams &params)
{
    params.numNodeItemIndices = 0;
    params.numWorkItems = 0;
    params.numActiveNodeItems = 0;
    numNodes = params.numNodes;
    numVMMs = params.numVMMs;
    ReserveBuffer(vmmUpdateNodeIndices, params.numVMMs);
    params.vmmUpdateNodeIndices = vmmUpdateNodeIndices.Data();
    UploadParams(params);
    Launch(*collectVMMUpdateNodes, 1, 1, "kd.collect_vmms");
    ReadBackParams(params);
}

void PathGuidingSystem::Update()
{
    ++renderedSamples;
    if (renderedSamples < nextUpdateSample)
        return;

    device->SyncDeviceStream();
    const auto updateStart = std::chrono::steady_clock::now();
    double lutMilliseconds = 0.0;
    stageTimers.clear();
    updateFrontierCount = 0;
    timingUpdate = true;
    DrainTrainingSamples();
    const size_t sampleCount = trainingSampleCount;
    if (sampleCount != 0)
    {
        YBI_ASSERT(sampleCount <= std::numeric_limits<uint32_t>::max());
        if (!vmfSplitTable.texture.textureObject)
        {
            const auto lutStart = std::chrono::steady_clock::now();
            vmfSplitTable = LoadVMFSplitTable(*device);
            lutMilliseconds = std::chrono::duration<double, std::milli>(
                                  std::chrono::steady_clock::now() - lutStart).count();
            if (!vmfSplitTable.texture.textureObject)
            {
                timingUpdate = false;
                return;
            }
        }
        PathGuidingKdTreeParams params{};
        vmfSplitTable.Bind(params);
        device->SyncDeviceStream();
        ReserveBuffer(trainingSampleIndices, sampleCount);
        ReserveBuffer(trainingSampleScratchIndices, sampleCount);
        params.sampleData = trainingSamples.Data();
        params.sampleIndices = trainingSampleIndices.Data();
        params.newSampleIndices = trainingSampleScratchIndices.Data();
        params.numSampleIndices = static_cast<uint32_t>(sampleCount);
        UploadParams(params);
        Launch(*initializeSampleIndices,
               CeilDiv(params.numSampleIndices, uint32_t(VMM_BLOCK_SIZE)),
               VMM_BLOCK_SIZE,
               "prepare.indices");
        UpdateKDTree(params);
        UpdateVMMs(params);
        device->SyncDeviceStream();
        trainingSampleCount = 0;
    }

    timingUpdate = false;
    const double wallMilliseconds = std::chrono::duration<double, std::milli>(
                                        std::chrono::steady_clock::now() - updateStart).count();
    ++iteration;
    PrintUpdateStatistics(sampleCount);
    PrintUpdateTimings(wallMilliseconds, lutMilliseconds);
    lastUpdateSample = renderedSamples;
    if (renderedSamples >= 8)
        updateSampleCount = Min(updateSampleCount * 2, 16u);
    nextUpdateSample = renderedSamples + updateSampleCount;
}

void PathGuidingSystem::PrintUpdateStatistics(size_t sampleCount) const
{
    const size_t fieldBytes = nodes.SizeInBytes() + vmmDistributions.SizeInBytes() +
                              vmmStatistics.SizeInBytes() + vmmSplitStatistics.SizeInBytes();
    const size_t totalBytes =
        fieldBytes + kernelParams.SizeInBytes() + radianceRecorderParams.SizeInBytes() +
        radianceRecords.SizeInBytes() +
        radianceRecordFreeList.SizeInBytes() + radianceRecordRetiredList.SizeInBytes() +
        radianceRecordAllocationCursorBuffer.SizeInBytes() +
        radianceRecordRetiredCountBuffer.SizeInBytes() +
        pathGuidingTrainingSamplesBuffer.SizeInBytes() +
        pathGuidingTrainingSampleCountBuffer.SizeInBytes() + trainingSamples.SizeInBytes() +
        trainingSampleIndices.SizeInBytes() + trainingSampleScratchIndices.SizeInBytes() +
        activeNodeItems.SizeInBytes() + nodeItems.SizeInBytes() + nodeItemIndices.SizeInBytes() +
        nextNodeItemIndices.SizeInBytes() + workItemLeftOffsets.SizeInBytes() +
        positionStatistics.SizeInBytes() + newPositionStatistics.SizeInBytes() +
        lookaheadStatistics.SizeInBytes() + newLookaheadStatistics.SizeInBytes() +
        guidingCellSplitMask.SizeInBytes() + vmmUpdateNodeIndices.SizeInBytes() +
        vmmSplitComponentMasks.SizeInBytes() +
        (vmfSplitTable.texture.textureObject ? VMFSplitTableSize * sizeof(float) : 0);
    constexpr double MiB = 1024.0 * 1024.0;
    std::printf("Path guiding iteration %u: spp=%llu, batch_spp=%llu, training_samples=%zu, "
                "dropped_total=%llu, nodes=%u, vmms=%u, field=%.2f MiB, total=%.2f MiB\n",
                iteration,
                static_cast<unsigned long long>(renderedSamples),
                static_cast<unsigned long long>(renderedSamples - lastUpdateSample),
                sampleCount,
                static_cast<unsigned long long>(droppedTrainingSampleCount),
                numNodes,
                numVMMs,
                fieldBytes / MiB,
                totalBytes / MiB);
    std::fflush(stdout);
}

void PathGuidingSystem::PrintUpdateTimings(double wallMilliseconds, double lutMilliseconds)
{
    struct StageTotal
    {
        double milliseconds = 0.0;
        uint32_t calls = 0;
    };
    std::map<std::string, StageTotal> totals;
    for (const StageTimer &entry : stageTimers)
    {
        StageTotal &total = totals[entry.stage];
        total.milliseconds += entry.timer->ElapsedMilliseconds();
        ++total.calls;
    }
    double kdMilliseconds = 0.0;
    double vmmMilliseconds = 0.0;
    double deviceMilliseconds = 0.0;
    for (const auto &entry : totals)
    {
        deviceMilliseconds += entry.second.milliseconds;
        if (entry.first.compare(0, 3, "kd.") == 0)
            kdMilliseconds += entry.second.milliseconds;
        if (entry.first.compare(0, 4, "vmm.") == 0)
            vmmMilliseconds += entry.second.milliseconds;
    }
    std::printf("Path guiding timing %u: spp=%llu, wall_ms=%.3f, device_ms=%.3f, "
                "kd_ms=%.3f, vmm_ms=%.3f, lut_ms=%.3f, frontiers=%u, timed_operations=%zu\n",
                iteration,
                static_cast<unsigned long long>(renderedSamples),
                wallMilliseconds, deviceMilliseconds, kdMilliseconds, vmmMilliseconds,
                lutMilliseconds, updateFrontierCount, stageTimers.size());
    for (const auto &entry : totals)
    {
        std::printf("Path guiding stage %u: %s, device_ms=%.3f, calls=%u\n",
                    iteration, entry.first.c_str(), entry.second.milliseconds, entry.second.calls);
    }
    stageTimers.clear();
    std::fflush(stdout);
}

void PathGuidingSystem::UpdateKDTree(PathGuidingKdTreeParams &params)
{
    BeginUpdate(params);
    while (params.numNodeItemIndices != 0)
    {
        if (timingUpdate)
            ++updateFrontierCount;
        ReserveFrontierBuffers(params);
        CalculatePositionStatistics(params);
        SelectNodeSplits(params);
        PartitionSamples(params);

        switch (splitDecision)
        {
            case PathGuidingSplitDecision::SampleCount:
                UpdateSampleCountSplits(params);
                break;
            case PathGuidingSplitDecision::Lookahead:
                UpdateLookaheadSplits(params);
                break;
        }
        AdvanceFrontier(params);
    }
    FinishUpdate(params);
}

void PathGuidingSystem::UpdateVMMs(PathGuidingKdTreeParams &params)
{
    if (params.numVMMUpdateNodes == 0)
        return;

    YBI_ASSERT(params.vmmUpdateNodeIndices && params.nodes && params.nodeItems &&
               params.vmmDistributions && params.vmmStatistics && params.vmmSplitStatistics);
    YBI_ASSERT(params.numSampleIndices == 0 || (params.sampleData && params.sampleIndices));
    ReserveVMMTrainingBuffers(params);
    UploadParams(params);
    Launch(*applyParallaxShift, params.numVMMUpdateNodes, 32, "vmm.parallax");
    Launch(*reprojectSamples, params.numVMMUpdateNodes, 32, "vmm.reproject");
    Launch(*updateVMM, params.numVMMUpdateNodes, VMM_BLOCK_SIZE, "vmm.fit");
    Launch(*splitVMM, params.numVMMUpdateNodes, VMM_BLOCK_SIZE, "vmm.split");
    Launch(*partialUpdateVMM, params.numVMMUpdateNodes, VMM_BLOCK_SIZE, "vmm.refit");
    Launch(*mergeVMM, params.numVMMUpdateNodes, VMM_BLOCK_SIZE, "vmm.merge");
    Launch(*updateComponentDistances, params.numVMMUpdateNodes, VMM_BLOCK_SIZE, "vmm.distance");
}

} // namespace ybi

#endif
