#include "kernel/backend/cuda/vmm_backend.h"
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

    CudaVMMBackend backend;
    if (backend.GetThreadIndex() == 0)
    {
        distribution.numComponents = params.distributions[mixtureIndex].numComponents;
        previousStats.numSamples = params.statistics[mixtureIndex].numSamples;
        YBI_ASSERT(distribution.numComponents > 0 &&
                   distribution.numComponents <= VMM_MAX_NUM_COMPONENTS);
    }
    backend.Sync();
    for (uint32_t k = backend.GetThreadIndex(); k < distribution.numComponents;
         k += backend.GetThreadStride())
    {
        distribution.weights[k] = params.distributions[mixtureIndex].weights[k];
        distribution.meanDirection[k] = params.distributions[mixtureIndex].meanDirection[k];
        distribution.kappa[k] = params.distributions[mixtureIndex].kappa[k];
        previousStats.sumWeightedDirections[k] =
            params.statistics[mixtureIndex].sumWeightedDirections[k];
        previousStats.sumWeightedStats[k] = params.statistics[mixtureIndex].sumWeightedStats[k];
    }
    backend.Sync();

    UpdateVMM(params.samples,
               distribution,
               previousStats,
               params.ranges[mixtureIndex],
               params.configuration,
               workspace,
               backend);

    for (uint32_t k = backend.GetThreadIndex(); k < distribution.numComponents;
         k += backend.GetThreadStride())
    {
        params.distributions[mixtureIndex].weights[k] = distribution.weights[k];
        params.distributions[mixtureIndex].meanDirection[k] = distribution.meanDirection[k];
        params.distributions[mixtureIndex].kappa[k] = distribution.kappa[k];
        params.statistics[mixtureIndex].sumWeightedDirections[k] =
            previousStats.sumWeightedDirections[k];
        params.statistics[mixtureIndex].sumWeightedStats[k] = previousStats.sumWeightedStats[k];
    }
    if (backend.GetThreadIndex() == 0)
    {
        params.statistics[mixtureIndex].numSamples = previousStats.numSamples;
    }
}

} // namespace ybi
