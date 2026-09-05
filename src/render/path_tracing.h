#pragma once

#include "core/frame.h"
#include "render/sampling.h"

namespace ybi
{

YBI_CPU_GPU float LinearToDisplay(float value)
{
    value = Clamp(value, 0.f, 1.f);
    return value <= 0.0031308f ? 12.92f * value : 1.055f * powf(value, 1.f / 2.4f) - 0.055f;
}

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

} // namespace ybi
