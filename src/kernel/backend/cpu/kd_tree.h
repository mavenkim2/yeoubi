#pragma once

#include "kernel/guiding/types.h"

namespace ybi
{

void CreateKdTreeNodesCPU(PathGuidingKdTreeParams &params, PathGuidingSplitDecision splitDecision);
void CreateKdTreeNodeItemsCPU(PathGuidingKdTreeParams &params, PathGuidingSplitDecision splitDecision);
void CreateKdTreeWorkItemsCPU(PathGuidingKdTreeParams &params, PathGuidingKdTreeStage stage);
void CalculateNodeBoundsCPU(PathGuidingKdTreeParams params);
void CalculateNodePositionStatisticsCPU(PathGuidingKdTreeParams params);
void CreateKdTreeReductionWorkItemsCPU(PathGuidingKdTreeParams &params);
void ReduceNodePositionStatisticsCPU(PathGuidingKdTreeParams params,
                                     const PathGuidingKdTreePositionStatistics *input,
                                     PathGuidingKdTreePositionStatistics *output);
void SelectNodeSplitsCPU(PathGuidingKdTreeParams params);
void PartitionCountCPU(PathGuidingKdTreeParams params);
void PartitionExclusiveScanCPU(PathGuidingKdTreeParams params);
void PartitionSamplesCPU(PathGuidingKdTreeParams params, uint32_t *newSampleIndices);

void CreateKdTreeNodesKernelCPU(void **arguments);
void CreateKdTreeNodeItemsKernelCPU(void **arguments);
void CreateKdTreeWorkItemsKernelCPU(void **arguments);
void CalculateNodeBoundsKernelCPU(void **arguments);
void CalculateNodePositionStatisticsKernelCPU(void **arguments);
void CreateKdTreeReductionWorkItemsKernelCPU(void **arguments);
void ReduceNodePositionStatisticsKernelCPU(void **arguments);
void SelectNodeSplitsKernelCPU(void **arguments);
void PartitionCountKernelCPU(void **arguments);
void PartitionExclusiveScanKernelCPU(void **arguments);
void PartitionSamplesKernelCPU(void **arguments);

} // namespace ybi
