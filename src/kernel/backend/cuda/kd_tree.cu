#include "kernel/backend/cuda/vmm_backend.h"
#include "kernel/backend/cuda/reduction.h"
#include "kernel/guiding/kd_tree.h"

namespace ybi
{

YBI_GPU void ValidateKdTreeBlockLaunch()
{
    YBI_ASSERT(blockDim.y == 1 && blockDim.z == 1);
    YBI_ASSERT(blockDim.x >= 32 && blockDim.x % 32 == 0);
    YBI_ASSERT(gridDim.y == 1 && gridDim.z == 1);
}

extern "C" __global__ void CreateKdTreeNodesKernel(PathGuidingKdTreeParams *kernelParams,
                                                  PathGuidingSplitDecision splitDecision)
{
    if (blockIdx.x == 0 && threadIdx.x == 0)
        CreateKdTreeNodes(*kernelParams, splitDecision);
}

extern "C" __global__ void CreateKdTreeNodeItemsKernel(PathGuidingKdTreeParams *kernelParams,
                                                      PathGuidingSplitDecision splitDecision)
{
    if (blockIdx.x == 0 && threadIdx.x == 0)
        CreateKdTreeNodeItems(*kernelParams, splitDecision);
}

extern "C" __global__ void CreateKdTreeWorkItemsKernel(PathGuidingKdTreeParams *kernelParams,
                                                      PathGuidingKdTreeStage stage)
{
    if (blockIdx.x == 0 && threadIdx.x == 0)
        CreateKdTreeWorkItems(*kernelParams, stage);
}

extern "C" __global__ void CalculateNodeBoundsKernel(const PathGuidingKdTreeParams *kernelParams)
{
    ValidateKdTreeBlockLaunch();
    const PathGuidingKdTreeParams &params = *kernelParams;
    for (uint32_t workItemIndex = blockIdx.x; workItemIndex < params.numWorkItems;
         workItemIndex += gridDim.x)
    {
        CalculateNodeBounds<CudaVMMBackend>(params, workItemIndex);
    }
}

extern "C" __global__ void CalculateNodePositionStatisticsKernel(
    const PathGuidingKdTreeParams *kernelParams)
{
    ValidateKdTreeBlockLaunch();
    const PathGuidingKdTreeParams &params = *kernelParams;
    for (uint32_t workItemIndex = blockIdx.x; workItemIndex < params.numWorkItems;
         workItemIndex += gridDim.x)
    {
        const PathGuidingKdTreeWorkItem workItem = GetKdTreeWorkItem(params, workItemIndex);
        const auto value = CalculateNodePositionStatistics(params, workItem, threadIdx.x);
        const auto statistics = BlockReduceSum<PathGuidingKdTreePositionStatistics, VMM_BLOCK_SIZE>(value);
        if (threadIdx.x == 0)
            StoreNodePositionStatistics(params, workItem.nodeItemIndex, workItemIndex,
                                        statistics, params.positionStatistics);
    }
}

extern "C" __global__ void CreateKdTreeReductionWorkItemsKernel(
    PathGuidingKdTreeParams *kernelParams)
{
    if (blockIdx.x == 0 && threadIdx.x == 0)
        CreateKdTreeReductionWorkItems(*kernelParams);
}

extern "C" __global__ void ReduceNodePositionStatisticsKernel(
    const PathGuidingKdTreeParams *kernelParams,
    const PathGuidingKdTreePositionStatistics *input,
    PathGuidingKdTreePositionStatistics *output)
{
    const PathGuidingKdTreeParams &params = *kernelParams;
    for (uint32_t workItemIndex = blockIdx.x; workItemIndex < params.numWorkItems;
         workItemIndex += gridDim.x)
    {
        const uint32_t nodeItemIndex = FindKdTreeNodeItemIndex(params, workItemIndex);
        const PathGuidingKdTreeNodeItem &nodeItem = params.nodeItems[nodeItemIndex];
        const auto value = LoadKdTreeReductionValue(nodeItem, workItemIndex, threadIdx.x, input);
        const auto statistics = BlockReduceSum<PathGuidingKdTreePositionStatistics, VMM_BLOCK_SIZE>(value);
        if (threadIdx.x == 0)
            StoreNodePositionStatistics(params, nodeItemIndex, workItemIndex, statistics, output);
    }
}

extern "C" __global__ void SelectNodeSplitsKernel(const PathGuidingKdTreeParams *kernelParams)
{
    YBI_ASSERT(blockDim.y == 1 && blockDim.z == 1);
    YBI_ASSERT(gridDim.y == 1 && gridDim.z == 1);
    const PathGuidingKdTreeParams &params = *kernelParams;
    const uint32_t threadIndex = blockIdx.x * blockDim.x + threadIdx.x;
    const uint32_t threadStride = blockDim.x * gridDim.x;
    for (uint32_t activeNodeItemIndex = threadIndex; activeNodeItemIndex < params.numActiveNodeItems;
         activeNodeItemIndex += threadStride)
    {
        SelectNodeSplits(params, params.activeNodeItems[activeNodeItemIndex]);
    }
}

extern "C" __global__ void PartitionCountKernel(const PathGuidingKdTreeParams *kernelParams)
{
    ValidateKdTreeBlockLaunch();
    const PathGuidingKdTreeParams &params = *kernelParams;
    for (uint32_t workItemIndex = blockIdx.x; workItemIndex < params.numWorkItems;
         workItemIndex += gridDim.x)
    {
        PartitionCount<CudaVMMBackend>(params, workItemIndex);
    }
}

extern "C" __global__ void PartitionExclusiveScanKernel(const PathGuidingKdTreeParams *kernelParams)
{
    ValidateKdTreeBlockLaunch();
    const PathGuidingKdTreeParams &params = *kernelParams;
    for (uint32_t activeNodeItemIndex = blockIdx.x; activeNodeItemIndex < params.numActiveNodeItems;
         activeNodeItemIndex += gridDim.x)
    {
        PartitionExclusiveScan<CudaVMMBackend>(params, params.activeNodeItems[activeNodeItemIndex]);
    }
}

extern "C" __global__ void PartitionSamplesKernel(const PathGuidingKdTreeParams *kernelParams,
                                                  uint32_t *newSampleIndices)
{
    ValidateKdTreeBlockLaunch();
    const PathGuidingKdTreeParams &params = *kernelParams;
    for (uint32_t workItemIndex = blockIdx.x; workItemIndex < params.numWorkItems;
         workItemIndex += gridDim.x)
    {
        PartitionSamples<CudaVMMBackend>(params, workItemIndex, newSampleIndices);
    }
}

} // namespace ybi
