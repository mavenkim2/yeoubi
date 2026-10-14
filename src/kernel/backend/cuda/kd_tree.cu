#include "kernel/backend/cuda/vmm_backend.h"
#include "kernel/guiding/kd_tree.h"

namespace ybi
{

YBI_GPU void ValidateKdTreeBlockLaunch()
{
    YBI_ASSERT(blockDim.y == 1 && blockDim.z == 1);
    YBI_ASSERT(blockDim.x >= 32 && blockDim.x % 32 == 0);
    YBI_ASSERT(gridDim.y == 1 && gridDim.z == 1);
}

extern "C" __global__ void CreateKdTreeNodeItemsKernel(const PathGuidingKdTreeParams *kernelParams)
{
    if (blockIdx.x == 0 && threadIdx.x == 0)
        CreateKdTreeNodeItems(*kernelParams);
}

extern "C" __global__ void CreateKdTreeWorkItemsKernel(const PathGuidingKdTreeParams *kernelParams)
{
    if (blockIdx.x == 0 && threadIdx.x == 0)
        CreateKdTreeWorkItems(*kernelParams);
}

extern "C" __global__ void CalculateNodeBoundsKernel(const PathGuidingKdTreeParams *kernelParams)
{
    ValidateKdTreeBlockLaunch();
    const PathGuidingKdTreeParams &params = *kernelParams;
    for (uint32_t workItemIndex = blockIdx.x; workItemIndex < params.numWorkItems;
         workItemIndex += gridDim.x)
    {
        CalculateNodeBounds<CudaVMMBackend>(params, workItemIndex);
    }
}

extern "C" __global__ void CalculateNodePositionStatisticsKernel(
    const PathGuidingKdTreeParams *kernelParams)
{
    ValidateKdTreeBlockLaunch();
    const PathGuidingKdTreeParams &params = *kernelParams;
    for (uint32_t workItemIndex = blockIdx.x; workItemIndex < params.numWorkItems;
         workItemIndex += gridDim.x)
    {
        CalculateNodePositionStatistics<CudaVMMBackend>(params, workItemIndex);
    }
}

extern "C" __global__ void PartitionCountKernel(const PathGuidingKdTreeParams *kernelParams)
{
    ValidateKdTreeBlockLaunch();
    const PathGuidingKdTreeParams &params = *kernelParams;
    for (uint32_t workItemIndex = blockIdx.x; workItemIndex < params.numWorkItems;
         workItemIndex += gridDim.x)
    {
        PartitionCount<CudaVMMBackend>(params, workItemIndex);
    }
}

extern "C" __global__ void PartitionExclusiveScanKernel(const PathGuidingKdTreeParams *kernelParams)
{
    YBI_ASSERT(blockDim.y == 1 && blockDim.z == 1);
    YBI_ASSERT(gridDim.y == 1 && gridDim.z == 1);
    const PathGuidingKdTreeParams &params = *kernelParams;
    const uint32_t threadIndex = blockIdx.x * blockDim.x + threadIdx.x;
    const uint32_t threadStride = blockDim.x * gridDim.x;
    for (uint32_t nodeItemIndex = threadIndex; nodeItemIndex < params.numNodeItems;
         nodeItemIndex += threadStride)
    {
        PartitionExclusiveScan(params, nodeItemIndex);
    }
}

extern "C" __global__ void PartitionSamplesKernel(const PathGuidingKdTreeParams *kernelParams,
                                                  uint32_t *newSampleIndices)
{
    ValidateKdTreeBlockLaunch();
    const PathGuidingKdTreeParams &params = *kernelParams;
    for (uint32_t workItemIndex = blockIdx.x; workItemIndex < params.numWorkItems;
         workItemIndex += gridDim.x)
    {
        PartitionSamples<CudaVMMBackend>(params, workItemIndex, newSampleIndices);
    }
}

} // namespace ybi
