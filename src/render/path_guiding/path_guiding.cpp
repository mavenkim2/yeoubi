#if defined(YBI_USE_PATH_GUIDING)

#if !defined(YBI_KD_TREE_KERNEL_PATH)
#error "YBI_USE_PATH_GUIDING requires YBI_KD_TREE_KERNEL_PATH"
#endif

#include "render/path_guiding/path_guiding.h"
#include "core/math.h"
#include "kernel/guiding/kd_tree.h"
#include <utility>

namespace ybi
{

PathGuidingSystem::PathGuidingSystem(Device *device, DeviceKind) : device(device)
{
    YBI_ASSERT(device);
    DeviceKernelDesc desc;
    desc.filename = YBI_KD_TREE_KERNEL_PATH;
    desc.entryName = "CreateKdTreeNodeItemsKernel";
    generateNodeItems = device->CreateDeviceKernel(desc);
    desc.entryName = "CreateKdTreeWorkItemsKernel";
    generateWorkItems = device->CreateDeviceKernel(desc);
    desc.entryName = "InitializeNodePositionStatisticsKernel";
    initializeNodePositionStatistics = device->CreateDeviceKernel(desc);
    desc.entryName = "CalculateNodeBoundsKernel";
    calculateNodeBounds = device->CreateDeviceKernel(desc);
    desc.entryName = "CalculateNodePositionStatisticsKernel";
    calculateNodePositionStatistics = device->CreateDeviceKernel(desc);
    desc.entryName = "SelectNodeSplitsKernel";
    selectNodeSplits = device->CreateDeviceKernel(desc);
    desc.entryName = "PartitionCountKernel";
    partitionCount = device->CreateDeviceKernel(desc);
    desc.entryName = "PartitionExclusiveScanKernel";
    partitionExclusiveScan = device->CreateDeviceKernel(desc);
    desc.entryName = "PartitionSamplesKernel";
    partitionSamples = device->CreateDeviceKernel(desc);
    YBI_ASSERT(generateNodeItems && generateWorkItems && initializeNodePositionStatistics &&
               calculateNodeBounds && calculateNodePositionStatistics && selectNodeSplits &&
               partitionCount && partitionExclusiveScan && partitionSamples);
    kernelParams = device->Alloc(sizeof(PathGuidingKdTreeParams));
    YBI_ASSERT(kernelParams);
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
    if (nextNodeItems)
        device->Free(nextNodeItems);
    if (workItems)
        device->Free(workItems);
}

void PathGuidingSystem::ReserveBuffer(DeviceBuffer &buffer, size_t size, size_t preserveSize)
{
    if (buffer.SizeInBytes() >= size)
        return;
    DeviceBuffer replacement = device->Alloc(Max(size, buffer.SizeInBytes() * 2));
    if (preserveSize != 0)
    {
        device->CopyDeviceToDeviceAsync(replacement.Data(), buffer.Data(), preserveSize);
        device->SyncDeviceStream();
    }
    device->Free(buffer);
    buffer = replacement;
}

void PathGuidingSystem::CopyToOwnedBuffer(DeviceBuffer &buffer, void *source, size_t size)
{
    const bool owned = buffer.Data() == source;
    ReserveBuffer(buffer, size, owned ? size : 0);
    if (!owned)
        device->CopyDeviceToDeviceAsync(buffer.Data(), source, size);
}

bool PathGuidingSystem::PartitionSamples(PathGuidingKdTreeParams &params)
{
    PathGuidingKdTreeStage stage = PathGuidingKdTreeStage::Partition;
    void *workItemArguments[] = {kernelParams.KernelArg(), &stage};
    DeviceKernelLaunchDesc desc;
    desc.arguments = workItemArguments;
    device->ExecuteKernel(*generateWorkItems, desc);
    device->SyncDeviceStream();
    device->CopyToHost(&params, kernelParams.Data(), sizeof(params));
    if (params.numWorkItems == 0)
        return false;

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
    return true;
}

void PathGuidingSystem::UpdateKDTree(PathGuidingKdTreeParams &params)
{
    params.generation = ++generation;
    params.numWorkItems = 0;
    params.numActiveNodeItems = 0;
    params.numNextNodeItems = 0;
    if (params.numNodeItems == 0)
        return;

    YBI_ASSERT(params.nodes && params.nodeItems && params.numNodes != 0);
    YBI_ASSERT(params.numSampleIndices == 0 || (params.sampleData && params.sampleIndices));
    CopyToOwnedBuffer(
        nodes, params.nodes, size_t(params.numNodes) * sizeof(PathGuidingKdTreeNode));
    CopyToOwnedBuffer(nodeItems,
                      params.nodeItems,
                      size_t(params.numNodeItems) * sizeof(PathGuidingKdTreeNodeItem));
    params.nodes = nodes.As<PathGuidingKdTreeNode>().Get();
    params.nodeItems = nodeItems.As<PathGuidingKdTreeNodeItem>().Get();

    void *arguments[] = {kernelParams.KernelArg()};
    PathGuidingKdTreeStage stage = PathGuidingKdTreeStage::PositionStatistics;
    void *workItemArguments[] = {kernelParams.KernelArg(), &stage};
    while (params.numNodeItems != 0)
    {
        ReserveBuffer(activeNodeItems, size_t(params.numNodeItems) * sizeof(uint32_t));
        const size_t workItemCapacity =
            CeilDiv(params.numSampleIndices, VMM_KD_TREE_WORK_ITEM_SAMPLE_LIMIT) +
            size_t(params.numNodeItems);
        ReserveBuffer(workItems, workItemCapacity * sizeof(PathGuidingKdTreeWorkItem));
        params.activeNodeItems = activeNodeItems.As<uint32_t>().Get();
        params.workItems = workItems.As<PathGuidingKdTreeWorkItem>().Get();
        stage = PathGuidingKdTreeStage::PositionStatistics;

        device->CopyToDevice(kernelParams.Data(), &params, sizeof(params));
        DeviceKernelLaunchDesc desc;
        desc.arguments = workItemArguments;
        device->ExecuteKernel(*generateWorkItems, desc);
        device->SyncDeviceStream();
        device->CopyToHost(&params, kernelParams.Data(), sizeof(params));
        if (params.numWorkItems != 0)
        {
            desc.arguments = arguments;
            desc.gridX = CeilDiv(params.numActiveNodeItems, VMM_BLOCK_SIZE);
            desc.blockX = VMM_BLOCK_SIZE;
            device->ExecuteKernel(*initializeNodePositionStatistics, desc);

            desc.gridX = params.numWorkItems;
            device->ExecuteKernel(*calculateNodeBounds, desc);
            device->ExecuteKernel(*calculateNodePositionStatistics, desc);

            // Select node splits.
            desc.gridX = CeilDiv(params.numActiveNodeItems, VMM_BLOCK_SIZE);
            device->ExecuteKernel(*selectNodeSplits, desc);
        }

        if (!PartitionSamples(params))
            break;

        ReserveBuffer(nodes,
                      (size_t(params.numNodes) + 2 * params.numActiveNodeItems) *
                          sizeof(PathGuidingKdTreeNode),
                      size_t(params.numNodes) * sizeof(PathGuidingKdTreeNode));
        ReserveBuffer(nextNodeItems,
                      size_t(2) * params.numActiveNodeItems * sizeof(PathGuidingKdTreeNodeItem));
        params.nodes = nodes.As<PathGuidingKdTreeNode>().Get();
        params.nextNodeItems = nextNodeItems.As<PathGuidingKdTreeNodeItem>().Get();
        device->CopyToDevice(kernelParams.Data(), &params, sizeof(params));
        DeviceKernelLaunchDesc nodeItemDesc;
        nodeItemDesc.arguments = arguments;
        device->ExecuteKernel(*generateNodeItems, nodeItemDesc);

        device->SyncDeviceStream();
        device->CopyToHost(&params, kernelParams.Data(), sizeof(params));
        std::swap(nodeItems, nextNodeItems);
        params.nodeItems = nodeItems.As<PathGuidingKdTreeNodeItem>().Get();
        params.nextNodeItems = nextNodeItems.As<PathGuidingKdTreeNodeItem>().Get();
        params.numNodeItems = params.numNextNodeItems;
        params.numNextNodeItems = 0;
    }
    params.numNodeItems = 0;
    params.numWorkItems = 0;
    params.numActiveNodeItems = 0;
    device->CopyToDevice(kernelParams.Data(), &params, sizeof(params));
}

} // namespace ybi

#endif
