#pragma once

#include "kernel/guiding/types.h"
#include "kernel/guiding/radiance_recorder.h"

namespace ybi
{

void ReturnRetiredRecordsKernel(RadianceRecorder *recorder);

/************************************************************
 * Kd-tree
 ************************************************************/

void InitializePathGuidingSampleIndicesKernel(PathGuidingKdTreeParams *kernelParams);
void CollectVMMUpdateNodesKernel(PathGuidingKdTreeParams *kernelParams);
void CreateKdTreeNodesKernel(PathGuidingKdTreeParams *kernelParams);
void CreateKdTreeNodeItemsKernel(PathGuidingKdTreeParams *kernelParams);
void CopySampleCountDistributionsKernel(PathGuidingKdTreeParams *kernelParams);
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
void PartialUpdateVMMKernel(PathGuidingKdTreeParams *kernelParams);
void SplitVMMKernel(PathGuidingKdTreeParams *kernelParams);
void MergeVMMKernel(PathGuidingKdTreeParams *kernelParams);

/************************************************************
 * Parallax
 ************************************************************/

void ReprojectSamplesKernel(PathGuidingKdTreeParams *kernelParams);
void UpdateComponentDistancesKernel(PathGuidingKdTreeParams *kernelParams);
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
