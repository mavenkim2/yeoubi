#include "kernel/state/types.h"
#include <cassert>
#include <cuda_runtime_api.h>

namespace ybi
{

extern "C" __global__ void StreamCompaction(WavefrontSampleState inputState,
                                            WavefrontSampleState outputState,
                                            uint32_t numInputSampleStates,
                                            uint32_t *numOutputSampleStates)
{
    assert(blockDim.x == 256 && blockDim.y == 1 && blockDim.z == 1);

    __shared__ uint32_t sharedMem[8];
    __shared__ uint32_t blockOffset;

    const uint32_t sampleIndex = blockDim.x * blockIdx.x + threadIdx.x;

    const uint32_t laneIndex = threadIdx.x & 31u;
    const uint32_t warpIndex = threadIdx.x >> 5u;

    const bool predicate = sampleIndex < numInputSampleStates &&
                           inputState.kernelTypes[sampleIndex] != KERNEL_FINISHED;

    const unsigned matches = __ballot_sync(0xffffffffu, predicate);
    const uint32_t count = __popc(matches);

    if (laneIndex == 0)
    {
        sharedMem[warpIndex] = count;
    }

    __syncthreads();

    if (threadIdx.x == 0)
    {
        uint32_t count = 0;
        for (uint32_t i = 0; i < 8; i++)
        {
            uint32_t value = sharedMem[i];
            sharedMem[i] = count;
            count += value;
        }

        blockOffset = count > 0 ? atomicAdd(numOutputSampleStates, count) : 0;
    }

    __syncthreads();

    const uint32_t rank = __popc(matches & ((1u << laneIndex) - 1u));
    const uint32_t laneOffset = blockOffset + sharedMem[warpIndex] + rank;

    if (predicate)
    {
        outputState.pixelIndex[laneOffset] = inputState.pixelIndex[sampleIndex];
        outputState.rng[laneOffset] = inputState.rng[sampleIndex];
        outputState.rayCone[laneOffset] = inputState.rayCone[sampleIndex];
        outputState.radiance[laneOffset] = inputState.radiance[sampleIndex];
#if defined(YBI_ENABLE_PATH_GUIDING)
        outputState.radianceRecordIndex[laneOffset] = inputState.radianceRecordIndex[sampleIndex];
        outputState.russianRouletteCorrectionFactor[laneOffset] =
            inputState.russianRouletteCorrectionFactor[sampleIndex];
#endif
        outputState.throughput[laneOffset] = inputState.throughput[sampleIndex];
        outputState.previousSurfacePosition[laneOffset] =
            inputState.previousSurfacePosition[sampleIndex];
        outputState.previousBsdfPdf[laneOffset] = inputState.previousBsdfPdf[sampleIndex];
        outputState.previousBsdfFlags[laneOffset] = inputState.previousBsdfFlags[sampleIndex];
        outputState.etaScale[laneOffset] = inputState.etaScale[sampleIndex];
        outputState.rayOrigin[laneOffset] = inputState.rayOrigin[sampleIndex];
        outputState.rayDirection[laneOffset] = inputState.rayDirection[sampleIndex];
        outputState.depth[laneOffset] = inputState.depth[sampleIndex];
        outputState.kernelTypes[laneOffset] = inputState.kernelTypes[sampleIndex];

        outputState.t[laneOffset] = inputState.t[sampleIndex];
        outputState.bary[laneOffset] = inputState.bary[sampleIndex];
        outputState.objectToWorld[laneOffset] = inputState.objectToWorld[sampleIndex];
        outputState.primIndex[laneOffset] = inputState.primIndex[sampleIndex];
        outputState.clusterId[laneOffset] = inputState.clusterId[sampleIndex];
        outputState.objectId[laneOffset] = inputState.objectId[sampleIndex];
        outputState.type[laneOffset] = inputState.type[sampleIndex];

        outputState.materialIndex[laneOffset] = inputState.materialIndex[sampleIndex];
        outputState.position[laneOffset] = inputState.position[sampleIndex];
        outputState.geometricNormal[laneOffset] = inputState.geometricNormal[sampleIndex];
        outputState.shadingNormal[laneOffset] = inputState.shadingNormal[sampleIndex];
        outputState.dpdu[laneOffset] = inputState.dpdu[sampleIndex];
        outputState.uv[laneOffset] = inputState.uv[sampleIndex];
        outputState.triangleArea[laneOffset] = inputState.triangleArea[sampleIndex];
        outputState.uvArea[laneOffset] = inputState.uvArea[sampleIndex];
        for (uint32_t i = 0; i < MAX_CLOSURE_SIZE; ++i)
        {
            outputState.closureData[size_t(laneOffset) * MAX_CLOSURE_SIZE + i] =
                inputState.closureData[size_t(sampleIndex) * MAX_CLOSURE_SIZE + i];
        }
    }
}

} // namespace ybi
