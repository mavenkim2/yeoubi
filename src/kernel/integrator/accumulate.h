#pragma once

#include "kernel/launch/kernel_params.h"
#include "kernel/integrator/path_tracing.h"

namespace ybi
{

YBI_GPU void Accumulate(const KernelParams &params, uint32_t pixelIndex)
{
    const Vec3f *kernelOutput = params.kernelOutput;
    Vec3f value = kernelOutput[pixelIndex];

    uint8_t *image = reinterpret_cast<uint8_t *>(params.imageHandle);
    Vec3f *history = params.accumulation ? params.accumulation + pixelIndex : nullptr;

    const Vec3f linearColor =
        AccumulateSamples(history, params.accumulatedSamples, value, params.samplesPerPixel);
    const Vec3f displayColor = FastACESToneMap(linearColor);

    image[4 * pixelIndex + 0] =
        static_cast<unsigned char>(255.f * LinearToDisplay(displayColor.x) + 0.5f);
    image[4 * pixelIndex + 1] =
        static_cast<unsigned char>(255.f * LinearToDisplay(displayColor.y) + 0.5f);
    image[4 * pixelIndex + 2] =
        static_cast<unsigned char>(255.f * LinearToDisplay(displayColor.z) + 0.5f);
    image[4 * pixelIndex + 3] = 255;
}

} // namespace ybi
