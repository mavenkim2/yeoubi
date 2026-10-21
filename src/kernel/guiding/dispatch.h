#pragma once

#include "kernel/guiding/kd_tree.h"
#include "kernel/guiding/lookahead.h"
#include "kernel/guiding/parallax.h"

#if !defined(YBI_KERNEL) || !defined(YBI_ITERATE_WORK_ITEMS) || !defined(YBI_ITERATE_ELEMENTS) || \
    !defined(YBI_END_ITERATE_WORK_ITEMS) || !defined(YBI_FIRST_THREAD) || \
    !defined(YBI_GUIDING_BACKEND)
#error "Guiding dispatch requires backend kernel macros and YBI_GUIDING_BACKEND"
#endif

namespace ybi
{

/************************************************************
 * Kd-tree
 ************************************************************/

YBI_KERNEL CreateKdTreeNodesKernel(PathGuidingKdTreeParams *kernelParams)
{
    if (YBI_FIRST_THREAD)
    {
        CreateKdTreeNodes(*kernelParams, kernelParams->splitDecision);
    }
}

YBI_KERNEL CreateKdTreeNodeItemsKernel(PathGuidingKdTreeParams *kernelParams)
{
    if (YBI_FIRST_THREAD)
    {
        CreateKdTreeNodeItems(*kernelParams, kernelParams->splitDecision);
    }
}

YBI_KERNEL CreateKdTreeWorkItemsKernel(PathGuidingKdTreeParams *kernelParams)
{
    if (YBI_FIRST_THREAD)
    {
        CreateKdTreeWorkItems(*kernelParams, kernelParams->stage);
    }
}

YBI_KERNEL CalculateNodeBoundsKernel(PathGuidingKdTreeParams *kernelParams)
{
    const PathGuidingKdTreeParams &params = *kernelParams;
    YBI_ITERATE_WORK_ITEMS(workItemIndex, params.numWorkItems)
    {
        CalculateNodeBounds<YBI_GUIDING_BACKEND>(params, workItemIndex);
    }
    YBI_END_ITERATE_WORK_ITEMS;
}

YBI_KERNEL CalculateNodePositionStatisticsKernel(PathGuidingKdTreeParams *kernelParams)
{
    const PathGuidingKdTreeParams &params = *kernelParams;
    YBI_ITERATE_WORK_ITEMS(workItemIndex, params.numWorkItems)
    {
        const PathGuidingKdTreeWorkItem workItem = GetKdTreeWorkItem(params, workItemIndex);
        const PathGuidingKdTreePositionStatistics statistics =
            YBI_GUIDING_BACKEND::ReduceBlock<PathGuidingKdTreePositionStatistics>(
                [&](uint32_t threadIndex) {
                    return CalculateNodePositionStatistics(params, workItem, threadIndex);
                });
        if (YBI_GUIDING_BACKEND::GetThreadIndex() == 0)
            params.positionStatistics[workItemIndex] = statistics;
    }
    YBI_END_ITERATE_WORK_ITEMS;
}

YBI_KERNEL CreateKdTreeReductionWorkItemsKernel(PathGuidingKdTreeParams *kernelParams)
{
    if (YBI_FIRST_THREAD)
    {
        CreateKdTreeReductionWorkItems(*kernelParams);
    }
}

YBI_KERNEL ReduceNodePositionStatisticsKernel(PathGuidingKdTreeParams *kernelParams)
{
    const PathGuidingKdTreeParams &params = *kernelParams;
    YBI_ITERATE_WORK_ITEMS(workItemIndex, params.numWorkItems)
    {
        const uint32_t nodeItemIndex = FindKdTreeNodeItemIndex(params, workItemIndex);
        const PathGuidingKdTreeNodeItem &nodeItem = params.nodeItems[nodeItemIndex];
        const PathGuidingKdTreePositionStatistics statistics =
            YBI_GUIDING_BACKEND::ReduceBlock<PathGuidingKdTreePositionStatistics>(
                [&](uint32_t threadIndex) {
                    return LoadKdTreeReductionValue(
                        nodeItem, workItemIndex, threadIndex, params.newPositionStatistics);
                });
        if (YBI_GUIDING_BACKEND::GetThreadIndex() == 0)
            params.positionStatistics[workItemIndex] = statistics;
    }
    YBI_END_ITERATE_WORK_ITEMS;
}

YBI_KERNEL AccumulateNodePositionStatisticsKernel(PathGuidingKdTreeParams *kernelParams)
{
    const PathGuidingKdTreeParams &params = *kernelParams;
    YBI_ITERATE_ELEMENTS(activeNodeItemIndex, params.numActiveNodeItems)
    {
        AccumulateNodePositionStatistics(
            params, params.activeNodeItems[activeNodeItemIndex], params.positionStatistics);
    }
    YBI_END_ITERATE_WORK_ITEMS;
}

YBI_KERNEL SelectNodeSplitsKernel(PathGuidingKdTreeParams *kernelParams)
{
    const PathGuidingKdTreeParams &params = *kernelParams;
    YBI_ITERATE_ELEMENTS(activeNodeItemIndex, params.numActiveNodeItems)
    {
        SelectNodeSplits(params, params.activeNodeItems[activeNodeItemIndex]);
    }
    YBI_END_ITERATE_WORK_ITEMS;
}

YBI_KERNEL PartitionCountKernel(PathGuidingKdTreeParams *kernelParams)
{
    const PathGuidingKdTreeParams &params = *kernelParams;
    YBI_ITERATE_WORK_ITEMS(workItemIndex, params.numWorkItems)
    {
        PartitionCount<YBI_GUIDING_BACKEND>(params, workItemIndex);
    }
    YBI_END_ITERATE_WORK_ITEMS;
}

YBI_KERNEL PartitionExclusiveScanKernel(PathGuidingKdTreeParams *kernelParams)
{
    const PathGuidingKdTreeParams &params = *kernelParams;
    YBI_ITERATE_WORK_ITEMS(activeNodeItemIndex, params.numActiveNodeItems)
    {
        PartitionExclusiveScan<YBI_GUIDING_BACKEND>(params,
                                                    params.activeNodeItems[activeNodeItemIndex]);
    }
    YBI_END_ITERATE_WORK_ITEMS;
}

YBI_KERNEL PartitionSamplesKernel(PathGuidingKdTreeParams *kernelParams)
{
    const PathGuidingKdTreeParams &params = *kernelParams;
    YBI_ITERATE_WORK_ITEMS(workItemIndex, params.numWorkItems)
    {
        PartitionSamples<YBI_GUIDING_BACKEND>(params, workItemIndex, params.newSampleIndices);
    }
    YBI_END_ITERATE_WORK_ITEMS;
}

/************************************************************
 * Parallax
 ************************************************************/

YBI_KERNEL ApplyParallaxShiftDispatchKernel(PathGuidingKdTreeParams *kernelParams)
{
    const PathGuidingKdTreeParams &params = *kernelParams;
    YBI_ITERATE_WORK_ITEMS(workItemIndex, params.numVMMUpdateNodes)
    {
        ApplyVMMParallaxShift<YBI_GUIDING_BACKEND>(params, workItemIndex);
    }
    YBI_END_ITERATE_WORK_ITEMS;
}

/************************************************************
 * Lookahead
 ************************************************************/

YBI_KERNEL UpdateLookaheadCellKernel(PathGuidingKdTreeParams *kernelParams)
{
    const PathGuidingKdTreeParams &params = *kernelParams;
    YBI_ITERATE_WORK_ITEMS(workItemIndex, params.numWorkItems)
    {
        const PathGuidingKdTreeWorkItem workItem = GetKdTreeWorkItem(params, workItemIndex);
        const PathGuidingKdTreeLookaheadStatistics statistics =
            YBI_GUIDING_BACKEND::ReduceBlock<PathGuidingKdTreeLookaheadStatistics>(
                [&](uint32_t threadIndex) {
                    return UpdateLookaheadCell(params, workItem, threadIndex);
                });
        if (YBI_GUIDING_BACKEND::GetThreadIndex() == 0)
            params.lookaheadStatistics[workItemIndex] = statistics;
    }
    YBI_END_ITERATE_WORK_ITEMS;
}

YBI_KERNEL ReduceLookaheadStatisticsKernel(PathGuidingKdTreeParams *kernelParams)
{
    const PathGuidingKdTreeParams &params = *kernelParams;
    YBI_ITERATE_WORK_ITEMS(workItemIndex, params.numWorkItems)
    {
        const uint32_t nodeItemIndex = FindKdTreeNodeItemIndex(params, workItemIndex);
        const PathGuidingKdTreeNodeItem &nodeItem = params.nodeItems[nodeItemIndex];
        const PathGuidingKdTreeLookaheadStatistics statistics =
            YBI_GUIDING_BACKEND::ReduceBlock<PathGuidingKdTreeLookaheadStatistics>(
                [&](uint32_t threadIndex) {
                    return LoadKdTreeReductionValue(
                        nodeItem, workItemIndex, threadIndex, params.newLookaheadStatistics);
                });
        if (YBI_GUIDING_BACKEND::GetThreadIndex() == 0)
            params.lookaheadStatistics[workItemIndex] = statistics;
    }
    YBI_END_ITERATE_WORK_ITEMS;
}

YBI_KERNEL AccumulateLookaheadCellsKernel(PathGuidingKdTreeParams *kernelParams)
{
    const PathGuidingKdTreeParams &params = *kernelParams;
    YBI_ITERATE_ELEMENTS(activeNodeItemIndex, params.numActiveNodeItems)
    {
        AccumulateLookaheadCell(
            params, params.activeNodeItems[activeNodeItemIndex], params.lookaheadStatistics);
    }
    YBI_END_ITERATE_WORK_ITEMS;
}

YBI_KERNEL ShouldSplitGuidingCellKernel(PathGuidingKdTreeParams *kernelParams)
{
    PathGuidingKdTreeParams &params = *kernelParams;
    YBI_ITERATE_ELEMENTS(activeNodeItemIndex, params.numActiveNodeItems)
    {
        ShouldSplitGuidingCell(params, params.activeNodeItems[activeNodeItemIndex]);
    }
    YBI_END_ITERATE_WORK_ITEMS;
}

YBI_KERNEL PromoteGuidingCellsKernel(PathGuidingKdTreeParams *kernelParams)
{
    if (YBI_FIRST_THREAD)
    {
        PromoteGuidingCells(*kernelParams);
    }
}

} // namespace ybi
