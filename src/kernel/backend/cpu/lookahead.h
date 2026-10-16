#pragma once

#include "kernel/guiding/types.h"

namespace ybi
{

void UpdateLookaheadCellCPU(PathGuidingKdTreeParams params);
void UpdateLookaheadCellKernelCPU(void **arguments);
void ReduceLookaheadStatisticsCPU(PathGuidingKdTreeParams params,
                                  const PathGuidingKdTreeLookaheadStatistics *input,
                                  PathGuidingKdTreeLookaheadStatistics *output);
void ReduceLookaheadStatisticsKernelCPU(void **arguments);
void AccumulateLookaheadCellsCPU(PathGuidingKdTreeParams params,
                                 const PathGuidingKdTreeLookaheadStatistics *input);
void AccumulateLookaheadCellsKernelCPU(void **arguments);
void ShouldSplitGuidingCellCPU(PathGuidingKdTreeParams &params);
void ShouldSplitGuidingCellKernelCPU(void **arguments);
void PromoteGuidingCellsCPU(PathGuidingKdTreeParams &params);
void PromoteGuidingCellsKernelCPU(void **arguments);

} // namespace ybi
