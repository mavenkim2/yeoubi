#pragma once

#include "kernel/guiding/types.h"

namespace ybi
{

void CreateKdTreeNodeItemsCPU(PathGuidingKdTreeParams params);
void CreateKdTreeWorkItemsCPU(PathGuidingKdTreeParams &params, PathGuidingKdTreeStage stage);
void InitializeNodePositionStatisticsCPU(PathGuidingKdTreeParams params);
void CalculateNodeBoundsCPU(PathGuidingKdTreeParams params);
void CalculateNodePositionStatisticsCPU(PathGuidingKdTreeParams params);
void PartitionCountCPU(PathGuidingKdTreeParams params);
void PartitionExclusiveScanCPU(PathGuidingKdTreeParams params);
void PartitionSamplesCPU(PathGuidingKdTreeParams params, uint32_t *newSampleIndices);

void CreateKdTreeNodeItemsKernelCPU(void **arguments);
void CreateKdTreeWorkItemsKernelCPU(void **arguments);
void InitializeNodePositionStatisticsKernelCPU(void **arguments);
void CalculateNodeBoundsKernelCPU(void **arguments);
void CalculateNodePositionStatisticsKernelCPU(void **arguments);
void PartitionCountKernelCPU(void **arguments);
void PartitionExclusiveScanKernelCPU(void **arguments);
void PartitionSamplesKernelCPU(void **arguments);

} // namespace ybi
