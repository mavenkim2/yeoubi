#include "kernel/backend/cuda/reduction.h"
#include "kernel/guiding/lookahead.h"

namespace ybi
{

extern "C" __global__ void UpdateLookaheadCellKernel(const PathGuidingKdTreeParams *kernelParams)
{
    YBI_ASSERT(blockDim.y == 1 && blockDim.z == 1);
    YBI_ASSERT(blockDim.x >= 32 && blockDim.x % 32 == 0);
    YBI_ASSERT(gridDim.y == 1 && gridDim.z == 1);
    const PathGuidingKdTreeParams &params = *kernelParams;
    for (uint32_t workItemIndex = blockIdx.x; workItemIndex < params.numWorkItems;
         workItemIndex += gridDim.x)
    {
        const PathGuidingKdTreeWorkItem workItem = GetKdTreeWorkItem(params, workItemIndex);
        const auto value = UpdateLookaheadCell(params, workItem, threadIdx.x);
        const auto statistics = BlockReduceSum<PathGuidingKdTreeLookaheadStatistics, VMM_BLOCK_SIZE>(value);
        if (threadIdx.x == 0)
            params.lookaheadStatistics[workItemIndex] = statistics;
    }
}

extern "C" __global__ void ReduceLookaheadStatisticsKernel(
    const PathGuidingKdTreeParams *kernelParams,
    const PathGuidingKdTreeLookaheadStatistics *input,
    PathGuidingKdTreeLookaheadStatistics *output)
{
    const PathGuidingKdTreeParams &params = *kernelParams;
    for (uint32_t workItemIndex = blockIdx.x; workItemIndex < params.numWorkItems;
         workItemIndex += gridDim.x)
    {
        const uint32_t nodeItemIndex = FindKdTreeNodeItemIndex(params, workItemIndex);
        const PathGuidingKdTreeNodeItem &nodeItem = params.nodeItems[nodeItemIndex];
        const auto value = LoadKdTreeReductionValue(nodeItem, workItemIndex, threadIdx.x, input);
        const auto statistics = BlockReduceSum<PathGuidingKdTreeLookaheadStatistics, VMM_BLOCK_SIZE>(value);
        if (threadIdx.x == 0)
            output[workItemIndex] = statistics;
    }
}

extern "C" __global__ void AccumulateLookaheadCellsKernel(
    const PathGuidingKdTreeParams *kernelParams,
    const PathGuidingKdTreeLookaheadStatistics *input)
{
    const PathGuidingKdTreeParams &params = *kernelParams;
    for (uint32_t activeIndex = blockIdx.x * blockDim.x + threadIdx.x;
         activeIndex < params.numActiveNodeItems; activeIndex += blockDim.x * gridDim.x)
        AccumulateLookaheadCell(params, params.activeNodeItems[activeIndex], input);
}

extern "C" __global__ void ShouldSplitGuidingCellKernel(PathGuidingKdTreeParams *kernelParams)
{
    YBI_ASSERT(blockDim.y == 1 && blockDim.z == 1);
    YBI_ASSERT(gridDim.y == 1 && gridDim.z == 1);
    PathGuidingKdTreeParams &params = *kernelParams;
    for (uint32_t activeNodeItemIndex = blockIdx.x * blockDim.x + threadIdx.x;
         activeNodeItemIndex < params.numActiveNodeItems;
         activeNodeItemIndex += blockDim.x * gridDim.x)
    {
        ShouldSplitGuidingCell(params, params.activeNodeItems[activeNodeItemIndex]);
    }
}

extern "C" __global__ void PromoteGuidingCellsKernel(PathGuidingKdTreeParams *kernelParams)
{
    if (blockIdx.x == 0 && threadIdx.x == 0)
        PromoteGuidingCells(*kernelParams);
}

} // namespace ybi
