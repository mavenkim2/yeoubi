#include "kernel/backend/cpu/kd_tree.h"
#include "kernel/backend/cpu/vmm_backend.h"
#include "kernel/guiding/kd_tree.h"
#include <cstring>
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

static PathGuidingKdTreeParams &GetKdTreeParams(void **arguments)
{
    YBI_ASSERT(arguments && arguments[0]);
    PathGuidingKdTreeParams *params = nullptr;
    std::memcpy(&params, arguments[0], sizeof(params));
    YBI_ASSERT(params);
    return *params;
}

void CreateKdTreeNodeItemsKernelCPU(void **arguments)
{
    CreateKdTreeNodeItemsCPU(GetKdTreeParams(arguments));
}

void CreateKdTreeWorkItemsKernelCPU(void **arguments)
{
    PathGuidingKdTreeParams &params = GetKdTreeParams(arguments);
    YBI_ASSERT(arguments[1]);
    const auto stage = *static_cast<PathGuidingKdTreeStage *>(arguments[1]);
    CreateKdTreeWorkItemsCPU(params, stage);
}

void InitializeNodePositionStatisticsKernelCPU(void **arguments)
{
    InitializeNodePositionStatisticsCPU(GetKdTreeParams(arguments));
}

void CalculateNodeBoundsKernelCPU(void **arguments)
{
    CalculateNodeBoundsCPU(GetKdTreeParams(arguments));
}

void CalculateNodePositionStatisticsKernelCPU(void **arguments)
{
    CalculateNodePositionStatisticsCPU(GetKdTreeParams(arguments));
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
