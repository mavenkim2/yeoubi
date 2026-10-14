#if defined(YBI_USE_PATH_GUIDING)

#if !defined(YBI_KD_TREE_KERNEL_PATH)
#error "YBI_USE_PATH_GUIDING requires YBI_KD_TREE_KERNEL_PATH"
#endif

#include "render/path_guiding/path_guiding.h"

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
    YBI_ASSERT(generateWorkItems && initializeNodePositionStatistics && calculateNodeBounds &&
               calculateNodePositionStatistics);
    kernelParams = device->Alloc(sizeof(PathGuidingKdTreeParams));
    YBI_ASSERT(kernelParams);
}

PathGuidingSystem::~PathGuidingSystem()
{
    if (kernelParams)
        device->Free(kernelParams);
}

void PathGuidingSystem::UpdateKDTree(PathGuidingKdTreeParams &params)
{
    params.numWorkItems = 0;
    if (params.numNodeItems == 0)
        return;

    YBI_ASSERT(params.nodes && params.nodeItems && params.workItems);
    YBI_ASSERT(params.sampleData && params.sampleIndices);
    PathGuidingKdTreeStage stage = PathGuidingKdTreeStage::PositionStatistics;

    device->CopyToDevice(kernelParams.Data(), &params, sizeof(params));
    void *workItemArguments[] = {kernelParams.KernelArg(), &stage};
    DeviceKernelLaunchDesc desc;
    desc.arguments = workItemArguments;
    device->ExecuteKernel(*generateWorkItems, desc);
    device->SyncDeviceStream();
    device->CopyToHost(&params, kernelParams.Data(), sizeof(params));
    if (params.numWorkItems == 0)
        return;

    void *arguments[] = {kernelParams.KernelArg()};
    desc.arguments = arguments;
    desc.gridX = (params.numNodeItems + VMM_BLOCK_SIZE - 1) / VMM_BLOCK_SIZE;
    desc.blockX = VMM_BLOCK_SIZE;
    device->ExecuteKernel(*initializeNodePositionStatistics, desc);

    desc.gridX = params.numWorkItems;
    device->ExecuteKernel(*calculateNodeBounds, desc);
    device->ExecuteKernel(*calculateNodePositionStatistics, desc);
    device->SyncDeviceStream();
}

} // namespace ybi

#endif
