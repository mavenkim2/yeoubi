#if defined(YBI_USE_PATH_GUIDING)

#if !defined(YBI_KD_TREE_KERNEL_PATH)
#error "YBI_USE_PATH_GUIDING requires YBI_KD_TREE_KERNEL_PATH"
#endif

#if !defined(YBI_LOOKAHEAD_KERNEL_PATH)
#error "YBI_USE_PATH_GUIDING requires YBI_LOOKAHEAD_KERNEL_PATH"
#endif

#include "render/path_guiding/path_guiding.h"
#include "core/math.h"
#include "kernel/guiding/kd_tree.h"
#include "kernel/guiding/lookahead.h"
#include <utility>

namespace ybi
{

PathGuidingSystem::PathGuidingSystem(Device *device,
                                     DeviceKind,
                                     PathGuidingSplitDecision splitDecision)
    : device(device), splitDecision(splitDecision)
{
    YBI_ASSERT(device);
    DeviceKernelDesc desc;
    desc.filename = YBI_KD_TREE_KERNEL_PATH;
    desc.entryName = "CreateKdTreeWorkItemsKernel";
    generateWorkItems = device->CreateDeviceKernel(desc);
    desc.entryName = "CalculateNodeBoundsKernel";
    calculateNodeBounds = device->CreateDeviceKernel(desc);
    desc.entryName = "CalculateNodePositionStatisticsKernel";
    calculateNodePositionStatistics = device->CreateDeviceKernel(desc);
    desc.entryName = "CreateKdTreeReductionWorkItemsKernel";
    generateReductionWorkItems = device->CreateDeviceKernel(desc);
    desc.entryName = "ReduceNodePositionStatisticsKernel";
    reduceNodePositionStatistics = device->CreateDeviceKernel(desc);
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
    YBI_ASSERT(generateWorkItems && calculateNodeBounds &&
               calculateNodePositionStatistics && generateReductionWorkItems &&
               reduceNodePositionStatistics && selectNodeSplits && partitionCount &&
               partitionExclusiveScan && partitionSamples && generateNodes && generateNodeItems);

    switch (splitDecision)
    {
        case PathGuidingSplitDecision::SampleCount:
            break;
        case PathGuidingSplitDecision::Lookahead:
            desc.filename = YBI_LOOKAHEAD_KERNEL_PATH;
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
            YBI_ASSERT(updateLookaheadCell && reduceLookaheadStatistics && accumulateLookaheadCells &&
                       shouldSplitGuidingCell && promoteGuidingCells);
            break;
    }
    kernelParams = device->Alloc<PathGuidingKdTreeParams>(1);
    YBI_ASSERT(kernelParams);

    PathGuidingKdTreeNode root{};
    if (splitDecision == PathGuidingSplitDecision::Lookahead)
    {
        root.guidingCellIndex = 0;
        VMMDistribution rootDistribution{};
        rootDistribution.numComponents = 1;
        rootDistribution.meanDirection[0] = Vec3f(0.f, 0.f, 1.f);
        rootDistribution.weights[0] = 1.f;
        vmmDistributions = device->Alloc<VMMDistribution>(1);
        device->CopyToDevice(vmmDistributions.Data(), &rootDistribution, sizeof(rootDistribution));
        numVMMs = 1;
    }
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

template <typename T>
T *PathGuidingSystem::ReduceNodeStatistics(const PathGuidingKdTreeParams &params, DeviceKernel &kernel,
                                           const DeviceBuffer<T> &input, const DeviceBuffer<T> &output)
{
    T *reductionInput = input.Data();
    T *reductionOutput = output.Data();
    void *arguments[] = {kernelParams.KernelArg()};
    void *reductionArguments[] = {kernelParams.KernelArg(), &reductionInput, &reductionOutput};
    uint32_t workItemBound = params.numWorkItems;
    DeviceKernelLaunchDesc desc;
    for (uint32_t count = params.maxWorkItemCount; count > 1; count = CeilDiv(count, VMM_BLOCK_SIZE))
    {
        desc.arguments = arguments;
        desc.gridX = 1;
        desc.blockX = 1;
        device->ExecuteKernel(*generateReductionWorkItems, desc);
        workItemBound = CeilDiv(workItemBound, VMM_BLOCK_SIZE) + params.numActiveNodeItems;
        desc.arguments = reductionArguments;
        desc.gridX = workItemBound;
        desc.blockX = VMM_BLOCK_SIZE;
        device->ExecuteKernel(kernel, desc);
        std::swap(reductionInput, reductionOutput);
    }
    return reductionInput;
}

void PathGuidingSystem::PartitionSamples(PathGuidingKdTreeParams &params)
{
    PathGuidingKdTreeStage stage = PathGuidingKdTreeStage::Partition;
    void *workItemArguments[] = {kernelParams.KernelArg(), &stage};
    DeviceKernelLaunchDesc desc;
    desc.arguments = workItemArguments;
    device->ExecuteKernel(*generateWorkItems, desc);
    device->SyncDeviceStream();
    device->CopyToHost(&params, kernelParams.Data(), sizeof(params));
    if (params.numWorkItems == 0)
        return;

    YBI_ASSERT(params.newSampleIndices && params.newSampleIndices != params.sampleIndices);
    void *arguments[] = {kernelParams.KernelArg()};
    desc.arguments = arguments;
    desc.blockX = VMM_BLOCK_SIZE;
    desc.gridX = params.numWorkItems;
    device->ExecuteKernel(*partitionCount, desc);

    desc.gridX = params.numActiveNodeItems;
    device->ExecuteKernel(*partitionExclusiveScan, desc);

    device->CopyDeviceToDeviceAsync(params.newSampleIndices,
                                    params.sampleIndices,
                                    size_t(params.numSampleIndices) * sizeof(uint32_t));
    void *partitionArguments[] = {kernelParams.KernelArg(), &params.newSampleIndices};
    desc.arguments = partitionArguments;
    desc.gridX = params.numWorkItems;
    device->ExecuteKernel(*partitionSamples, desc);
    std::swap(params.sampleIndices, params.newSampleIndices);
}

void PathGuidingSystem::CreateNodeItems(PathGuidingKdTreeParams &params)
{
    void *arguments[] = {kernelParams.KernelArg(), &splitDecision};
    const size_t nodeCapacity = size_t(params.numNodes) + size_t(2) * params.numNodeItemIndices;
    ReserveBuffer(nodes, nodeCapacity, params.numNodes);
    ReserveBuffer(nodeItems, nodeCapacity, params.numNodes);
    ReserveBuffer(nextNodeItemIndices,
                  size_t(2) * (size_t(params.numNodeItemIndices) + params.numGuidingCellSplits));
    params.nodes = nodes.Data();
    params.nodeItems = nodeItems.Data();
    params.nextNodeItemIndices = nextNodeItemIndices.Data();
    device->CopyToDevice(kernelParams.Data(), &params, sizeof(params));
    DeviceKernelLaunchDesc nodeItemDesc;
    nodeItemDesc.arguments = arguments;
    device->ExecuteKernel(*generateNodes, nodeItemDesc);
    device->ExecuteKernel(*generateNodeItems, nodeItemDesc);

    device->SyncDeviceStream();
    device->CopyToHost(&params, kernelParams.Data(), sizeof(params));
}

void PathGuidingSystem::UpdateSampleCountSplits(PathGuidingKdTreeParams &params)
{
    CreateNodeItems(params);
}

void PathGuidingSystem::UpdateLookaheadSplits(PathGuidingKdTreeParams &params)
{
    const size_t nodeCapacity = size_t(params.numNodes) + size_t(2) * params.numNodeItemIndices;
    const size_t maskWordCount = CeilDiv(nodeCapacity, size_t(32));
    ReserveBuffer(guidingCellSplitMask, maskWordCount);
    device->Memzero(guidingCellSplitMask.Data(), maskWordCount * sizeof(uint32_t));
    params.guidingCellSplitMask = guidingCellSplitMask.Data();
    params.numGuidingCellSplits = 0;
    device->CopyToDevice(kernelParams.Data(), &params, sizeof(params));

    PathGuidingKdTreeStage stage = PathGuidingKdTreeStage::Lookahead;
    void *workItemArguments[] = {kernelParams.KernelArg(), &stage};
    DeviceKernelLaunchDesc desc;
    desc.arguments = workItemArguments;
    device->ExecuteKernel(*generateWorkItems, desc);
    device->SyncDeviceStream();
    device->CopyToHost(&params, kernelParams.Data(), sizeof(params));
    void *arguments[] = {kernelParams.KernelArg()};
    desc.arguments = arguments;
    if (params.numWorkItems != 0)
    {
        ReserveBuffer(lookaheadStatistics, params.numWorkItems);
        ReserveBuffer(newLookaheadStatistics, params.numWorkItems);
        params.lookaheadStatistics = lookaheadStatistics.Data();
        device->CopyToDevice(kernelParams.Data(), &params, sizeof(params));
        desc.blockX = VMM_BLOCK_SIZE;
        desc.gridX = params.numWorkItems;
        device->ExecuteKernel(*updateLookaheadCell, desc);
        auto *statistics = ReduceNodeStatistics(params, *reduceLookaheadStatistics,
                                                 lookaheadStatistics, newLookaheadStatistics);
        void *accumulateArguments[] = {kernelParams.KernelArg(), &statistics};
        desc.arguments = accumulateArguments;
        desc.gridX = CeilDiv(params.numActiveNodeItems, VMM_BLOCK_SIZE);
        device->ExecuteKernel(*accumulateLookaheadCells, desc);
        desc.arguments = arguments;
        desc.gridX = CeilDiv(params.numActiveNodeItems, VMM_BLOCK_SIZE);
        device->ExecuteKernel(*shouldSplitGuidingCell, desc);
        device->SyncDeviceStream();
        device->CopyToHost(&params, kernelParams.Data(), sizeof(params));
    }

    CreateNodeItems(params);

    if (params.numGuidingCellSplits != 0)
    {
        ReserveBuffer(vmmDistributions,
                      size_t(params.numVMMs) + params.numGuidingCellSplits,
                      params.numVMMs);
        params.vmmDistributions = vmmDistributions.Data();
        device->CopyToDevice(kernelParams.Data(), &params, sizeof(params));
        desc.gridX = 1;
        desc.blockX = 1;
        device->ExecuteKernel(*promoteGuidingCells, desc);
        device->SyncDeviceStream();
        device->CopyToHost(&params, kernelParams.Data(), sizeof(params));
    }
}

void PathGuidingSystem::UpdateKDTree(PathGuidingKdTreeParams &params)
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
    params.numVMMs = numVMMs;

    void *arguments[] = {kernelParams.KernelArg()};
    PathGuidingKdTreeStage stage = PathGuidingKdTreeStage::PositionStatistics;
    void *workItemArguments[] = {kernelParams.KernelArg(), &stage};
    while (params.numNodeItemIndices != 0)
    {
        ReserveBuffer(activeNodeItems, params.numNodeItemIndices);
        const size_t workItemCapacity =
            CeilDiv(params.numSampleIndices, VMM_KD_TREE_WORK_ITEM_SAMPLE_LIMIT) +
            size_t(params.numNodeItemIndices);
        ReserveBuffer(workItemLeftOffsets, workItemCapacity);
        params.activeNodeItems = activeNodeItems.Data();
        params.workItemLeftOffsets = workItemLeftOffsets.Data();
        stage = PathGuidingKdTreeStage::PositionStatistics;

        device->CopyToDevice(kernelParams.Data(), &params, sizeof(params));
        DeviceKernelLaunchDesc desc;
        desc.arguments = workItemArguments;
        device->ExecuteKernel(*generateWorkItems, desc);
        device->SyncDeviceStream();
        device->CopyToHost(&params, kernelParams.Data(), sizeof(params));
        if (params.numWorkItems != 0)
        {
            ReserveBuffer(positionStatistics, params.numWorkItems);
            ReserveBuffer(newPositionStatistics, params.numWorkItems);
            params.positionStatistics = positionStatistics.Data();
            device->CopyToDevice(kernelParams.Data(), &params, sizeof(params));
            desc.arguments = arguments;
            desc.blockX = VMM_BLOCK_SIZE;
            desc.gridX = params.numWorkItems;
            device->ExecuteKernel(*calculateNodeBounds, desc);
            device->ExecuteKernel(*calculateNodePositionStatistics, desc);

            ReduceNodeStatistics(params, *reduceNodePositionStatistics,
                                 positionStatistics, newPositionStatistics);

            // Select node splits.
            desc.arguments = arguments;
            desc.gridX = CeilDiv(params.numActiveNodeItems, VMM_BLOCK_SIZE);
            device->ExecuteKernel(*selectNodeSplits, desc);
        }

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
        params.numNodeItemIndices = params.numNextNodeItemIndices;
        params.numNextNodeItemIndices = 0;
        std::swap(nodeItemIndices, nextNodeItemIndices);
        params.nodeItemIndices = nodeItemIndices.Data();
        params.nextNodeItemIndices = nextNodeItemIndices.Data();
    }
    params.numNodeItemIndices = 0;
    params.numWorkItems = 0;
    params.numActiveNodeItems = 0;
    numNodes = params.numNodes;
    numVMMs = params.numVMMs;
    device->CopyToDevice(kernelParams.Data(), &params, sizeof(params));
}

} // namespace ybi

#endif
