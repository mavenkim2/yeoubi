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

YBI_CPU_GPU float ReinhardToneMap(float value)
{
    value = Max(value, 0.f);
    return value / (1.f + value);
}

YBI_CPU_GPU Vec3f ReinhardToneMap(const Vec3f &value)
{
    return Vec3f(ReinhardToneMap(value.x), ReinhardToneMap(value.y), ReinhardToneMap(value.z));
}

YBI_CPU_GPU float FastACESToneMap(float value)
{
    value = Max(value, 0.f);

    constexpr float a = 2.51f;
    constexpr float b = 0.03f;
    constexpr float c = 2.43f;
    constexpr float d = 0.59f;
    constexpr float e = 0.14f;

    return Clamp((value * (a * value + b)) / (value * (c * value + d) + e), 0.f, 1.f);
}

YBI_CPU_GPU Vec3f FastACESToneMap(const Vec3f &value)
{
    return Vec3f(FastACESToneMap(value.x), FastACESToneMap(value.y), FastACESToneMap(value.z));
}

YBI_CPU_GPU float PowerHeuristic(float sampledPdf, float otherPdf)
{
    const float sampledPdfSquared = sampledPdf * sampledPdf;
    const float otherPdfSquared = otherPdf * otherPdf;
    const float sum = sampledPdfSquared + otherPdfSquared;
    return sum > 0.f ? sampledPdfSquared / sum : 0.f;
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
