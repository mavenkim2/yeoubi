#include "kernel/backend/cpu/lookahead.h"
#include "kernel/backend/cpu/reduction.h"
#include "kernel/guiding/lookahead.h"
#include <cstring>
#include <tbb/parallel_for.h>

namespace ybi
{

void UpdateLookaheadCellCPU(PathGuidingKdTreeParams params)
{
    tbb::parallel_for(0u, params.numWorkItems, [&](uint32_t workItemIndex) {
        const PathGuidingKdTreeWorkItem workItem = GetKdTreeWorkItem(params, workItemIndex);
        PathGuidingKdTreeLookaheadStatistics values[VMM_BLOCK_SIZE] = {};
        for (uint32_t thread = 0; thread < VMM_BLOCK_SIZE; ++thread)
            values[thread] = UpdateLookaheadCell(params, workItem, thread);
        params.lookaheadStatistics[workItemIndex] =
            BlockReduceSum<PathGuidingKdTreeLookaheadStatistics, VMM_BLOCK_SIZE>(values);
    });
}

void ReduceLookaheadStatisticsCPU(PathGuidingKdTreeParams params,
                                  const PathGuidingKdTreeLookaheadStatistics *input,
                                  PathGuidingKdTreeLookaheadStatistics *output)
{
    tbb::parallel_for(0u, params.numWorkItems, [&](uint32_t workItemIndex) {
        const uint32_t nodeItemIndex = FindKdTreeNodeItemIndex(params, workItemIndex);
        const PathGuidingKdTreeNodeItem &nodeItem = params.nodeItems[nodeItemIndex];
        PathGuidingKdTreeLookaheadStatistics values[VMM_BLOCK_SIZE] = {};
        for (uint32_t thread = 0; thread < VMM_BLOCK_SIZE; ++thread)
            values[thread] = LoadKdTreeReductionValue(nodeItem, workItemIndex, thread, input);
        output[workItemIndex] =
            BlockReduceSum<PathGuidingKdTreeLookaheadStatistics, VMM_BLOCK_SIZE>(values);
    });
}

void AccumulateLookaheadCellsCPU(PathGuidingKdTreeParams params,
                                 const PathGuidingKdTreeLookaheadStatistics *input)
{
    tbb::parallel_for(0u, params.numActiveNodeItems, [&](uint32_t activeNodeItemIndex) {
        AccumulateLookaheadCell(params, params.activeNodeItems[activeNodeItemIndex], input);
    });
}

void ShouldSplitGuidingCellCPU(PathGuidingKdTreeParams &params)
{
    tbb::parallel_for(0u, params.numActiveNodeItems, [&](uint32_t activeNodeItemIndex) {
        ShouldSplitGuidingCell(params, params.activeNodeItems[activeNodeItemIndex]);
    });
}

static PathGuidingKdTreeParams &GetLookaheadParams(void **arguments)
{
    YBI_ASSERT(arguments && arguments[0]);
    PathGuidingKdTreeParams *params = nullptr;
    std::memcpy(&params, arguments[0], sizeof(params));
    YBI_ASSERT(params);
    return *params;
}

void UpdateLookaheadCellKernelCPU(void **arguments)
{
    UpdateLookaheadCellCPU(GetLookaheadParams(arguments));
}

void ReduceLookaheadStatisticsKernelCPU(void **arguments)
{
    const PathGuidingKdTreeLookaheadStatistics *input = nullptr;
    PathGuidingKdTreeLookaheadStatistics *output = nullptr;
    std::memcpy(&input, arguments[1], sizeof(input));
    std::memcpy(&output, arguments[2], sizeof(output));
    ReduceLookaheadStatisticsCPU(GetLookaheadParams(arguments), input, output);
}

void AccumulateLookaheadCellsKernelCPU(void **arguments)
{
    const PathGuidingKdTreeLookaheadStatistics *input = nullptr;
    std::memcpy(&input, arguments[1], sizeof(input));
    AccumulateLookaheadCellsCPU(GetLookaheadParams(arguments), input);
}

void ShouldSplitGuidingCellKernelCPU(void **arguments)
{
    ShouldSplitGuidingCellCPU(GetLookaheadParams(arguments));
}

void PromoteGuidingCellsCPU(PathGuidingKdTreeParams &params)
{
    PromoteGuidingCells(params);
}

void PromoteGuidingCellsKernelCPU(void **arguments)
{
    PromoteGuidingCellsCPU(GetLookaheadParams(arguments));
}

} // namespace ybi
