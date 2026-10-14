#include "render/path_guiding/path_guiding.h"
#include "kernel/backend/cpu/kd_tree.h"
#include <stdexcept>

namespace ybi
{

PathGuidingSystem::PathGuidingSystem(Device *device, DeviceKind deviceKind)
    : device(device), deviceKind(deviceKind)
{
    YBI_ASSERT(device);
    if (deviceKind == DeviceKind::CPU)
        return;

#if defined(YBI_KD_TREE_KERNEL_PATH)
    generateWorkItems =
        device->CreateDeviceKernel({YBI_KD_TREE_KERNEL_PATH, "CreateKdTreeWorkItemsKernel"});
    initializeNodePositionStatistics = device->CreateDeviceKernel(
        {YBI_KD_TREE_KERNEL_PATH, "InitializeNodePositionStatisticsKernel"});
    calculateNodeBounds =
        device->CreateDeviceKernel({YBI_KD_TREE_KERNEL_PATH, "CalculateNodeBoundsKernel"});
    calculateNodePositionStatistics = device->CreateDeviceKernel(
        {YBI_KD_TREE_KERNEL_PATH, "CalculateNodePositionStatisticsKernel"});
    kernelParams = device->Alloc(sizeof(PathGuidingKdTreeParams));
#else
    throw std::runtime_error("CUDA kd-tree kernels are unavailable");
#endif
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

    if (deviceKind == DeviceKind::CPU)
    {
#if defined(SIMD_4)
        CreateKdTreeWorkItemsCPU(params, stage);
        if (params.numWorkItems == 0)
            return;
        InitializeNodePositionStatisticsCPU(params);
        CalculateNodeBoundsCPU(params);
        CalculateNodePositionStatisticsCPU(params);
        return;
#else
        throw std::runtime_error("CPU kd-tree helpers require SIMD_4");
#endif
    }

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
