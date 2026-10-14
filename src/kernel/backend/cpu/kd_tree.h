#pragma once

#include "kernel/guiding/types.h"

namespace ybi
{

void CreateKdTreeNodeItemsCPU(PathGuidingKdTreeParams params);
void CreateKdTreeWorkItemsCPU(PathGuidingKdTreeParams params);
void CalculateNodeBoundsCPU(PathGuidingKdTreeParams params);
void CalculateNodePositionStatisticsCPU(PathGuidingKdTreeParams params);
void PartitionCountCPU(PathGuidingKdTreeParams params);
void PartitionExclusiveScanCPU(PathGuidingKdTreeParams params);
void PartitionSamplesCPU(PathGuidingKdTreeParams params, uint32_t *newSampleIndices);

} // namespace ybi
