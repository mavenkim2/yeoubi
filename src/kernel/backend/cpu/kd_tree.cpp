#include "kernel/backend/cpu/kd_tree.h"
#include "kernel/backend/cpu/vmm_backend.h"
#include "kernel/guiding/kd_tree.h"
#include <tbb/parallel_for.h>

namespace ybi
{

void CreateKdTreeNodeItemsCPU(PathGuidingKdTreeParams params)
{
    CreateKdTreeNodeItems(params);
}

void CreateKdTreeWorkItemsCPU(PathGuidingKdTreeParams &params, PathGuidingKdTreeStage stage)
{
    CreateKdTreeWorkItems(params, stage);
}

void InitializeNodePositionStatisticsCPU(PathGuidingKdTreeParams params)
{
    tbb::parallel_for(0u, params.numNodeItems, [&](uint32_t nodeItemIndex) {
        InitializeNodePositionStatistics(params, nodeItemIndex);
    });
}

void CalculateNodeBoundsCPU(PathGuidingKdTreeParams params)
{
    tbb::parallel_for(0u, params.numNodeItems, [&](uint32_t nodeItemIndex) {
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
    tbb::parallel_for(0u, params.numNodeItems, [&](uint32_t nodeItemIndex) {
        const PathGuidingKdTreeNodeItem &nodeItem = params.nodeItems[nodeItemIndex];
        for (uint32_t workItemIndex = nodeItem.workItemOffset;
             workItemIndex < nodeItem.workItemOffset + nodeItem.workItemCount;
             ++workItemIndex)
        {
            CalculateNodePositionStatistics<CpuSIMDVMMBackend>(params, workItemIndex);
        }
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
    tbb::parallel_for(0u, params.numNodeItems, [&](uint32_t nodeItemIndex) {
        PartitionExclusiveScan(params, nodeItemIndex);
    });
}

void PartitionSamplesCPU(PathGuidingKdTreeParams params, uint32_t *newSampleIndices)
{
    tbb::parallel_for(0u, params.numWorkItems, [&](uint32_t workItemIndex) {
        PartitionSamples<CpuSIMDVMMBackend>(params, workItemIndex, newSampleIndices);
    });
}

} // namespace ybi
