#pragma once

#include "kernel/film/color.h"
#include "kernel/launch/kernel_params.h"

namespace ybi
{

YBI_CPU_GPU Vec3f AccumulateSamples(Vec3f *history,
                                    uint32_t previousSamples,
                                    const Vec3f &batchSum,
                                    uint32_t batchSamples)
{
    const Vec3f batchMean = batchSum / float(batchSamples);
    if (!history)
        return batchMean;
    // Reset frames overwrite history without reading uninitialized or stale pixels.
    if (previousSamples == 0)
        *history = batchMean;
    else
        *history +=
            (batchMean - *history) * (float(batchSamples) / float(previousSamples + batchSamples));
    return *history;
}

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
