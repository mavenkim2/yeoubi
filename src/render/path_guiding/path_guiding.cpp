#if defined(YBI_USE_PATH_GUIDING)

#if !defined(YBI_KD_TREE_KERNEL_PATH)
#error "YBI_USE_PATH_GUIDING requires YBI_KD_TREE_KERNEL_PATH"
#endif

#include "render/path_guiding/path_guiding.h"
#include <utility>

namespace ybi
{

PathGuidingSystem::PathGuidingSystem(Device *device, DeviceKind) : device(device)
{
    YBI_ASSERT(device);
    DeviceKernelDesc desc;
    desc.filename = YBI_KD_TREE_KERNEL_PATH;
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
    YBI_ASSERT(generateWorkItems && initializeNodePositionStatistics && calculateNodeBounds &&
               calculateNodePositionStatistics && selectNodeSplits && partitionCount &&
               partitionExclusiveScan && partitionSamples);
    kernelParams = device->Alloc(sizeof(PathGuidingKdTreeParams));
    YBI_ASSERT(kernelParams);
}

PathGuidingSystem::~PathGuidingSystem()
{
    if (kernelParams)
        device->Free(kernelParams);
    if (activeNodeItems)
        device->Free(activeNodeItems);
}

void PathGuidingSystem::UpdateKDTree(PathGuidingKdTreeParams &params)
{
    params.generation = ++generation;
    params.numWorkItems = 0;
    params.numActiveNodeItems = 0;
    if (params.numNodeItems == 0)
        return;

    YBI_ASSERT(params.nodes && params.nodeItems && params.workItems);
    YBI_ASSERT(params.sampleData && params.sampleIndices);
    const size_t activeNodeItemBytes = size_t(params.numNodeItems) * sizeof(uint32_t);
    if (activeNodeItems.SizeInBytes() < activeNodeItemBytes)
    {
        device->Free(activeNodeItems);
        activeNodeItems = device->Alloc(activeNodeItemBytes);
    }
    params.activeNodeItems = activeNodeItems.As<uint32_t>().Get();
    PathGuidingKdTreeStage stage = PathGuidingKdTreeStage::PositionStatistics;

    device->CopyToDevice(kernelParams.Data(), &params, sizeof(params));
    void *workItemArguments[] = {kernelParams.KernelArg(), &stage};
    DeviceKernelLaunchDesc desc;
    desc.arguments = workItemArguments;
    device->ExecuteKernel(*generateWorkItems, desc);
    device->SyncDeviceStream();
    device->CopyToHost(&params, kernelParams.Data(), sizeof(params));
    void *arguments[] = {kernelParams.KernelArg()};
    if (params.numWorkItems != 0)
    {
        desc.arguments = arguments;
        desc.gridX = (params.numActiveNodeItems + VMM_BLOCK_SIZE - 1) / VMM_BLOCK_SIZE;
        desc.blockX = VMM_BLOCK_SIZE;
        device->ExecuteKernel(*initializeNodePositionStatistics, desc);

        desc.gridX = params.numWorkItems;
        device->ExecuteKernel(*calculateNodeBounds, desc);
        device->ExecuteKernel(*calculateNodePositionStatistics, desc);

        // Select node splits.
        desc.gridX = (params.numActiveNodeItems + VMM_BLOCK_SIZE - 1) / VMM_BLOCK_SIZE;
        device->ExecuteKernel(*selectNodeSplits, desc);
    }

    // Partition samples.
    stage = PathGuidingKdTreeStage::Partition;
    desc = DeviceKernelLaunchDesc{};
    desc.arguments = workItemArguments;
    device->ExecuteKernel(*generateWorkItems, desc);
    device->SyncDeviceStream();
    device->CopyToHost(&params, kernelParams.Data(), sizeof(params));
    if (params.numWorkItems == 0)
        return;

    YBI_ASSERT(params.newSampleIndices && params.newSampleIndices != params.sampleIndices);
    desc.arguments = arguments;
    desc.blockX = VMM_BLOCK_SIZE;
    desc.gridX = params.numWorkItems;
    device->ExecuteKernel(*partitionCount, desc);

    desc.gridX = (params.numActiveNodeItems + VMM_BLOCK_SIZE - 1) / VMM_BLOCK_SIZE;
    device->ExecuteKernel(*partitionExclusiveScan, desc);

    device->CopyDeviceToDeviceAsync(params.newSampleIndices,
                                    params.sampleIndices,
                                    size_t(params.numSampleIndices) * sizeof(uint32_t));

    void *partitionArguments[] = {kernelParams.KernelArg(), &params.newSampleIndices};
    desc.arguments = partitionArguments;
    desc.gridX = params.numWorkItems;
    device->ExecuteKernel(*partitionSamples, desc);
    device->SyncDeviceStream();
    std::swap(params.sampleIndices, params.newSampleIndices);
    device->CopyToDevice(kernelParams.Data(), &params, sizeof(params));
}

} // namespace ybi

#endif
