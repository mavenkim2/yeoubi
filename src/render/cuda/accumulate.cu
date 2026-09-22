#include "device/kernel_params.h"
#include "render/accumulator.h"

namespace ybi
{

extern "C" __global__ void AccumulateSamples(KernelParams params)
{
    const uint32_t x = blockDim.x * blockIdx.x + threadIdx.x;
    const uint32_t y = blockDim.y * blockIdx.y + threadIdx.y;

    if (x >= static_cast<uint32_t>(params.imageWidth) ||
        y >= static_cast<uint32_t>(params.imageHeight))
        return;

    const uint32_t pixelIndex = x + y * static_cast<uint32_t>(params.imageWidth);
    Accumulate(params, pixelIndex);
}
} // namespace ybi
