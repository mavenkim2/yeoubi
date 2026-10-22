#pragma once

#include "kernel/guiding/types.h"

namespace ybi
{

/************************************************************
 * Kd-tree
 ************************************************************/

void CollectVMMUpdateNodesKernel(PathGuidingKdTreeParams *kernelParams);
void CreateKdTreeNodesKernel(PathGuidingKdTreeParams *kernelParams);
void CreateKdTreeNodeItemsKernel(PathGuidingKdTreeParams *kernelParams);
void CreateKdTreeWorkItemsKernel(PathGuidingKdTreeParams *kernelParams);
void CalculateNodeBoundsKernel(PathGuidingKdTreeParams *kernelParams);
void CalculateNodePositionStatisticsKernel(PathGuidingKdTreeParams *kernelParams);
void CreateKdTreeReductionWorkItemsKernel(PathGuidingKdTreeParams *kernelParams);
void ReduceNodePositionStatisticsKernel(PathGuidingKdTreeParams *kernelParams);
void AccumulateNodePositionStatisticsKernel(PathGuidingKdTreeParams *kernelParams);
void SelectNodeSplitsKernel(PathGuidingKdTreeParams *kernelParams);
void PartitionCountKernel(PathGuidingKdTreeParams *kernelParams);
void PartitionExclusiveScanKernel(PathGuidingKdTreeParams *kernelParams);
void PartitionSamplesKernel(PathGuidingKdTreeParams *kernelParams);

/************************************************************
 * VMM
 ************************************************************/

void UpdateVMMKernel(PathGuidingKdTreeParams *kernelParams);

/************************************************************
 * Parallax
 ************************************************************/

void ReprojectSamplesKernel(PathGuidingKdTreeParams *kernelParams);
void ApplyParallaxShiftDispatchKernel(PathGuidingKdTreeParams *kernelParams);

/************************************************************
 * Lookahead
 ************************************************************/

void UpdateLookaheadCellKernel(PathGuidingKdTreeParams *kernelParams);
void ReduceLookaheadStatisticsKernel(PathGuidingKdTreeParams *kernelParams);
void AccumulateLookaheadCellsKernel(PathGuidingKdTreeParams *kernelParams);
void ShouldSplitGuidingCellKernel(PathGuidingKdTreeParams *kernelParams);
void PromoteGuidingCellsKernel(PathGuidingKdTreeParams *kernelParams);

} // namespace ybi
