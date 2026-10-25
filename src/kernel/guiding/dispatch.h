#pragma once

#include "kernel/guiding/guiding.h"
#include "kernel/guiding/kd_tree.h"
#include "kernel/guiding/lookahead.h"
#include "kernel/guiding/parallax.h"
#include "kernel/guiding/radiance_recorder.h"

#if !defined(YBI_KERNEL) || !defined(YBI_SHARED) || !defined(YBI_ITERATE_WORK_ITEMS) || \
    !defined(YBI_ITERATE_ELEMENTS) || !defined(YBI_END_ITERATE_WORK_ITEMS) || \
    !defined(YBI_FIRST_THREAD) || !defined(YBI_GUIDING_BACKEND)
#error "Guiding dispatch requires backend kernel macros and YBI_GUIDING_BACKEND"
#endif

namespace ybi
{

YBI_KERNEL ReturnRetiredRecordsKernel(RadianceRecorder *recorder)
{
    YBI_ITERATE_WORK_ITEMS(workItemIndex, 1u)
    {
        ReturnRetiredRecords<YBI_GUIDING_BACKEND>(*recorder);
    }
    YBI_END_ITERATE_WORK_ITEMS;
}

/************************************************************
 * Kd-tree
 ************************************************************/

YBI_KERNEL InitializePathGuidingSampleIndicesKernel(PathGuidingKdTreeParams *kernelParams)
{
    const PathGuidingKdTreeParams &params = *kernelParams;
    YBI_ITERATE_ELEMENTS(sampleIndex, params.numSampleIndices)
    {
        params.sampleIndices[sampleIndex] = sampleIndex;
    }
    YBI_END_ITERATE_WORK_ITEMS;
}

YBI_KERNEL CollectVMMUpdateNodesKernel(PathGuidingKdTreeParams *kernelParams)
{
    if (YBI_FIRST_THREAD)
    {
        CollectVMMUpdateNodes(*kernelParams);
    }
}

YBI_KERNEL CreateKdTreeNodesKernel(PathGuidingKdTreeParams *kernelParams)
{
    if (YBI_FIRST_THREAD)
    {
        CreateKdTreeNodes(*kernelParams);
    }
}

YBI_KERNEL CreateKdTreeNodeItemsKernel(PathGuidingKdTreeParams *kernelParams)
{
    if (YBI_FIRST_THREAD)
    {
        CreateKdTreeNodeItems(*kernelParams, kernelParams->splitDecision);
    }
}

YBI_KERNEL CopySampleCountDistributionsKernel(PathGuidingKdTreeParams *kernelParams)
{
    if (YBI_FIRST_THREAD)
    {
        CopySampleCountDistributions(*kernelParams);
    }
}

YBI_KERNEL CreateKdTreeWorkItemsKernel(PathGuidingKdTreeParams *kernelParams)
{
    YBI_ITERATE_WORK_ITEMS(workItemIndex, 1u)
    {
        CreateKdTreeWorkItems<YBI_GUIDING_BACKEND>(*kernelParams, kernelParams->stage);
    }
    YBI_END_ITERATE_WORK_ITEMS;
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
 * VMM
 ************************************************************/

YBI_KERNEL UpdateVMMKernel(PathGuidingKdTreeParams *kernelParams)
{
    const PathGuidingKdTreeParams &params = *kernelParams;
    YBI_ITERATE_WORK_ITEMS(updateIndex, params.numVMMUpdateNodes)
    {
        YBI_SHARED VMMDistribution distribution;
        YBI_SHARED SufficientStatistics previousStats;
        YBI_SHARED VMMUpdateWorkspace workspace;

        UpdateNodeVMM<YBI_GUIDING_BACKEND>(
            params, updateIndex, distribution, previousStats, workspace);
    }
    YBI_END_ITERATE_WORK_ITEMS;
}

YBI_KERNEL PartialUpdateVMMKernel(PathGuidingKdTreeParams *kernelParams)
{
    const PathGuidingKdTreeParams &params = *kernelParams;
    YBI_ITERATE_WORK_ITEMS(updateIndex, params.numVMMUpdateNodes)
    {
        YBI_SHARED VMMDistribution distribution;
        YBI_SHARED SufficientStatistics previousStats;
        YBI_SHARED VMMUpdateWorkspace workspace;

        UpdateNodeVMM<YBI_GUIDING_BACKEND, true>(
            params, updateIndex, distribution, previousStats, workspace);
    }
    YBI_END_ITERATE_WORK_ITEMS;
}

YBI_KERNEL SplitVMMKernel(PathGuidingKdTreeParams *kernelParams)
{
    const PathGuidingKdTreeParams &params = *kernelParams;
    YBI_ITERATE_WORK_ITEMS(updateIndex, params.numVMMUpdateNodes)
    {
        YBI_SHARED VMMDistribution distribution;
        YBI_SHARED SufficientStatistics previousStats;
        YBI_SHARED VMMSplitWorkspace workspace;

        SplitNodeVMM<YBI_GUIDING_BACKEND>(
            params, updateIndex, distribution, previousStats, workspace);
    }
    YBI_END_ITERATE_WORK_ITEMS;
}

YBI_KERNEL MergeVMMKernel(PathGuidingKdTreeParams *kernelParams)
{
    const PathGuidingKdTreeParams &params = *kernelParams;
    YBI_ITERATE_WORK_ITEMS(updateIndex, params.numVMMUpdateNodes)
    {
        YBI_SHARED VMMDistribution distribution;
        YBI_SHARED SufficientStatistics previousStats;
        YBI_SHARED VMMMergeWorkspace workspace;

        MergeNodeVMM<YBI_GUIDING_BACKEND>(
            params, updateIndex, distribution, previousStats, workspace);
    }
    YBI_END_ITERATE_WORK_ITEMS;
}

/************************************************************
 * Parallax
 ************************************************************/

YBI_KERNEL UpdateComponentDistancesKernel(PathGuidingKdTreeParams *kernelParams)
{
    const PathGuidingKdTreeParams &params = *kernelParams;
    YBI_ITERATE_WORK_ITEMS(updateIndex, params.numVMMUpdateNodes)
    {
        YBI_SHARED VMMDistribution distribution;
        YBI_SHARED SufficientStatistics previousStats;
        YBI_SHARED VMMComponentDistanceWorkspace workspace;

        UpdateNodeComponentDistances<YBI_GUIDING_BACKEND>(
            params, updateIndex, distribution, previousStats, workspace);
    }
    YBI_END_ITERATE_WORK_ITEMS;
}

YBI_KERNEL ReprojectSamplesKernel(PathGuidingKdTreeParams *kernelParams)
{
    const PathGuidingKdTreeParams &params = *kernelParams;
    YBI_ITERATE_WORK_ITEMS(workItemIndex, params.numVMMUpdateNodes)
    {
        ReprojectSamples<YBI_GUIDING_BACKEND>(params, workItemIndex);
    }
    YBI_END_ITERATE_WORK_ITEMS;
}

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
    YBI_ITERATE_WORK_ITEMS(workItemIndex, 1u)
    {
        YBI_SHARED uint32_t vmmOffset;
        YBI_SHARED uint32_t frontierOffset;
        PromoteGuidingCells<YBI_GUIDING_BACKEND>(
            *kernelParams, vmmOffset, frontierOffset);
    }
    YBI_END_ITERATE_WORK_ITEMS;
}

} // namespace ybi
