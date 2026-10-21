#include "kernel/backend/cuda/guiding/vmm_backend.h"
#include "kernel/launch/vmm_update_params.h"

namespace ybi
{

// Launch one block per mixture with a one-dimensional, warp-aligned block
// (e.g. 256 threads). Each mixture's range must contain at least one sample.
extern "C" __global__ void UpdateVMMKernel(const VMMUpdateParams *kernelParams)
{
    YBI_ASSERT(blockDim.y == 1 && blockDim.z == 1);
    YBI_ASSERT(blockDim.x >= 32 && blockDim.x % 32 == 0);
    YBI_ASSERT(gridDim.y == 1 && gridDim.z == 1);

    const VMMUpdateParams &params = *kernelParams;
    const uint32_t mixtureIndex = blockIdx.x;
    if (mixtureIndex >= params.numDistributions)
    {
        return;
    }

    __shared__ VMMDistribution distribution;
    __shared__ SufficientStatistics previousStats;
    __shared__ VMMUpdateWorkspace workspace;

    using Backend = CudaVMMBackend;
    if (Backend::GetThreadIndex() == 0)
    {
        distribution.numComponents = params.distributions[mixtureIndex].numComponents;
        distribution.pivotPosition = params.distributions[mixtureIndex].pivotPosition;
        distribution.initialized = params.distributions[mixtureIndex].initialized;
        previousStats.numSamples = params.statistics[mixtureIndex].numSamples;
        YBI_ASSERT(distribution.numComponents > 0 &&
                   distribution.numComponents <= VMM_MAX_NUM_COMPONENTS);
    }
    Backend::Sync();
    for (uint32_t k = Backend::GetThreadIndex(); k < distribution.numComponents;
         k += Backend::GetThreadStride())
    {
        distribution.weights[k] = params.distributions[mixtureIndex].weights[k];
        distribution.meanDirection[k] = params.distributions[mixtureIndex].meanDirection[k];
        distribution.kappa[k] = params.distributions[mixtureIndex].kappa[k];
        distribution.distance[k] = params.distributions[mixtureIndex].distance[k];
        previousStats.sumWeightedDirections[k] =
            params.statistics[mixtureIndex].sumWeightedDirections[k];
        previousStats.sumWeightedStats[k] = params.statistics[mixtureIndex].sumWeightedStats[k];
    }
    Backend::Sync();

    UpdateVMM<Backend>(params.samples,
                       distribution,
                       previousStats,
                       params.ranges[mixtureIndex],
                       params.configuration,
                       workspace);

    for (uint32_t k = Backend::GetThreadIndex(); k < distribution.numComponents;
         k += Backend::GetThreadStride())
    {
        params.distributions[mixtureIndex].weights[k] = distribution.weights[k];
        params.distributions[mixtureIndex].meanDirection[k] = distribution.meanDirection[k];
        params.distributions[mixtureIndex].kappa[k] = distribution.kappa[k];
        params.distributions[mixtureIndex].distance[k] = distribution.distance[k];
        params.statistics[mixtureIndex].sumWeightedDirections[k] =
            previousStats.sumWeightedDirections[k];
        params.statistics[mixtureIndex].sumWeightedStats[k] = previousStats.sumWeightedStats[k];
    }
    if (Backend::GetThreadIndex() == 0)
    {
        params.distributions[mixtureIndex].pivotPosition = distribution.pivotPosition;
        params.distributions[mixtureIndex].initialized = distribution.initialized;
        params.statistics[mixtureIndex].numSamples = previousStats.numSamples;
    }
}

} // namespace ybi
