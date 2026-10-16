#include "kernel/backend/cpu/kd_tree.h"
#include "kernel/backend/cpu/reduction.h"
#include "kernel/backend/cpu/vmm_backend.h"
#include "kernel/guiding/kd_tree.h"
#include <cstring>
#include <tbb/parallel_for.h>

namespace ybi
{

void CreateKdTreeNodesCPU(PathGuidingKdTreeParams &params, PathGuidingSplitDecision splitDecision)
{
    CreateKdTreeNodes(params, splitDecision);
}

void CreateKdTreeNodeItemsCPU(PathGuidingKdTreeParams &params, PathGuidingSplitDecision splitDecision)
{
    CreateKdTreeNodeItems(params, splitDecision);
}

void CreateKdTreeWorkItemsCPU(PathGuidingKdTreeParams &params, PathGuidingKdTreeStage stage)
{
    CreateKdTreeWorkItems(params, stage);
}

void CalculateNodeBoundsCPU(PathGuidingKdTreeParams params)
{
    tbb::parallel_for(0u, params.numActiveNodeItems, [&](uint32_t activeNodeItemIndex) {
        const uint32_t nodeItemIndex = params.activeNodeItems[activeNodeItemIndex];
        const PathGuidingKdTreeNodeItem &nodeItem = params.nodeItems[nodeItemIndex];
        for (uint32_t workItemIndex = nodeItem.workItemOffset;
             workItemIndex < nodeItem.workItemOffset + nodeItem.workItemCount;
             ++workItemIndex)
        {
            CalculateNodeBounds<CpuSIMDVMMBackend>(params, workItemIndex);
        }
    });
}

void CalculateNodePositionStatisticsCPU(PathGuidingKdTreeParams params)
{
    tbb::parallel_for(0u, params.numWorkItems, [&](uint32_t workItemIndex) {
        const PathGuidingKdTreeWorkItem workItem = GetKdTreeWorkItem(params, workItemIndex);
        PathGuidingKdTreePositionStatistics values[VMM_BLOCK_SIZE] = {};
        for (uint32_t thread = 0; thread < VMM_BLOCK_SIZE; ++thread)
            values[thread] = CalculateNodePositionStatistics(params, workItem, thread);
        const auto statistics = BlockReduceSum<PathGuidingKdTreePositionStatistics, VMM_BLOCK_SIZE>(values);
        StoreNodePositionStatistics(params, workItem.nodeItemIndex, workItemIndex,
                                    statistics, params.positionStatistics);
    });
}

void CreateKdTreeReductionWorkItemsCPU(PathGuidingKdTreeParams &params)
{
    CreateKdTreeReductionWorkItems(params);
}

void ReduceNodePositionStatisticsCPU(PathGuidingKdTreeParams params,
                                     const PathGuidingKdTreePositionStatistics *input,
                                     PathGuidingKdTreePositionStatistics *output)
{
    tbb::parallel_for(0u, params.numWorkItems, [&](uint32_t workItemIndex) {
        const uint32_t nodeItemIndex = FindKdTreeNodeItemIndex(params, workItemIndex);
        const PathGuidingKdTreeNodeItem &nodeItem = params.nodeItems[nodeItemIndex];
        PathGuidingKdTreePositionStatistics values[VMM_BLOCK_SIZE] = {};
        for (uint32_t thread = 0; thread < VMM_BLOCK_SIZE; ++thread)
            values[thread] = LoadKdTreeReductionValue(nodeItem, workItemIndex, thread, input);
        const auto statistics = BlockReduceSum<PathGuidingKdTreePositionStatistics, VMM_BLOCK_SIZE>(values);
        StoreNodePositionStatistics(params, nodeItemIndex, workItemIndex, statistics, output);
    });
}

void SelectNodeSplitsCPU(PathGuidingKdTreeParams params)
{
    tbb::parallel_for(0u, params.numActiveNodeItems, [&](uint32_t activeNodeItemIndex) {
        SelectNodeSplits(params, params.activeNodeItems[activeNodeItemIndex]);
    });
}

void PartitionCountCPU(PathGuidingKdTreeParams params)
{
    tbb::parallel_for(0u, params.numWorkItems, [&](uint32_t workItemIndex) {
        PartitionCount<CpuSIMDVMMBackend>(params, workItemIndex);
    });
}

void PartitionExclusiveScanCPU(PathGuidingKdTreeParams params)
{
    tbb::parallel_for(0u, params.numActiveNodeItems, [&](uint32_t activeNodeItemIndex) {
        PartitionExclusiveScan<CpuSIMDVMMBackend>(params, params.activeNodeItems[activeNodeItemIndex]);
    });
}

void PartitionSamplesCPU(PathGuidingKdTreeParams params, uint32_t *newSampleIndices)
{
    tbb::parallel_for(0u, params.numWorkItems, [&](uint32_t workItemIndex) {
        PartitionSamples<CpuSIMDVMMBackend>(params, workItemIndex, newSampleIndices);
    });
}

static PathGuidingKdTreeParams &GetKdTreeParams(void **arguments)
{
    YBI_ASSERT(arguments && arguments[0]);
    PathGuidingKdTreeParams *params = nullptr;
    std::memcpy(&params, arguments[0], sizeof(params));
    YBI_ASSERT(params);
    return *params;
}

void CreateKdTreeNodesKernelCPU(void **arguments)
{
    YBI_ASSERT(arguments && arguments[1]);
    const auto splitDecision = *static_cast<PathGuidingSplitDecision *>(arguments[1]);
    CreateKdTreeNodesCPU(GetKdTreeParams(arguments), splitDecision);
}

void CreateKdTreeNodeItemsKernelCPU(void **arguments)
{
    YBI_ASSERT(arguments && arguments[1]);
    const auto splitDecision = *static_cast<PathGuidingSplitDecision *>(arguments[1]);
    CreateKdTreeNodeItemsCPU(GetKdTreeParams(arguments), splitDecision);
}

void CreateKdTreeWorkItemsKernelCPU(void **arguments)
{
    PathGuidingKdTreeParams &params = GetKdTreeParams(arguments);
    YBI_ASSERT(arguments[1]);
    const auto stage = *static_cast<PathGuidingKdTreeStage *>(arguments[1]);
    CreateKdTreeWorkItemsCPU(params, stage);
}

void CalculateNodeBoundsKernelCPU(void **arguments)
{
    CalculateNodeBoundsCPU(GetKdTreeParams(arguments));
}

void CalculateNodePositionStatisticsKernelCPU(void **arguments)
{
    CalculateNodePositionStatisticsCPU(GetKdTreeParams(arguments));
}

void CreateKdTreeReductionWorkItemsKernelCPU(void **arguments)
{
    CreateKdTreeReductionWorkItemsCPU(GetKdTreeParams(arguments));
}

void ReduceNodePositionStatisticsKernelCPU(void **arguments)
{
    const PathGuidingKdTreePositionStatistics *input = nullptr;
    PathGuidingKdTreePositionStatistics *output = nullptr;
    std::memcpy(&input, arguments[1], sizeof(input));
    std::memcpy(&output, arguments[2], sizeof(output));
    ReduceNodePositionStatisticsCPU(GetKdTreeParams(arguments), input, output);
}

void SelectNodeSplitsKernelCPU(void **arguments)
{
    SelectNodeSplitsCPU(GetKdTreeParams(arguments));
}

void PartitionCountKernelCPU(void **arguments)
{
    PartitionCountCPU(GetKdTreeParams(arguments));
}

void PartitionExclusiveScanKernelCPU(void **arguments)
{
    PartitionExclusiveScanCPU(GetKdTreeParams(arguments));
}

void PartitionSamplesKernelCPU(void **arguments)
{
    PathGuidingKdTreeParams &params = GetKdTreeParams(arguments);
    YBI_ASSERT(arguments[1]);
    uint32_t *newSampleIndices = nullptr;
    std::memcpy(&newSampleIndices, arguments[1], sizeof(newSampleIndices));
    YBI_ASSERT(newSampleIndices);
    PartitionSamplesCPU(params, newSampleIndices);
}

} // namespace ybi
