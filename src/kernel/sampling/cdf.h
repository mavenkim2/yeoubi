#pragma once

#include "core/math.h"
#include "core/vec2f.h"
#include "device/cuda_compat.h"

namespace ybi
{

YBI_CPU_GPU bool SamplePiecewise1DCDF(const float *cdf,
                                      uint32_t numValues,
                                      float u,
                                      uint32_t &index,
                                      float &fraction,
                                      float &probability)
{
    u = Clamp(u, 0.f, FLOAT_BEFORE_ONE);
    uint32_t low = 0;
    uint32_t high = numValues;
    while (low + 1 < high)
    {
        const uint32_t middle = low + (high - low) / 2;
        if (cdf[middle] <= u)
            low = middle;
        else
            high = middle;
    }

    index = low;
    probability = cdf[index + 1] - cdf[index];
    YBI_ASSERT(probability != 0.f);

    fraction = Clamp((u - cdf[index]) / probability, 0.f, FLOAT_BEFORE_ONE);
    return true;
}

YBI_CPU_GPU bool SamplePiecewise2DCDF(const float *conditionalCDF,
                                      const float *marginalCDF,
                                      uint32_t width,
                                      uint32_t height,
                                      const Vec2f &u,
                                      Vec2f &uv,
                                      float &pdf)
{
    pdf = 0.f;
    uint32_t x, y;
    float dx, dy;
    float columnProbability, rowProbability;
    if (!SamplePiecewise1DCDF(marginalCDF, height, u.y, y, dy, rowProbability))
        return false;
    const float *row = conditionalCDF + size_t(y) * (size_t(width) + 1);
    if (!SamplePiecewise1DCDF(row, width, u.x, x, dx, columnProbability))
        return false;

    uv = Vec2f(Min((float(x) + dx) / float(width), FLOAT_BEFORE_ONE),
               Min((float(y) + dy) / float(height), FLOAT_BEFORE_ONE));
    pdf = rowProbability * columnProbability * float(width) * float(height);
    return pdf > 0.f;
}

YBI_CPU_GPU float PDFPiecewise2DCDF(const float *conditionalCDF,
                                    const float *marginalCDF,
                                    uint32_t width,
                                    uint32_t height,
                                    const Vec2f &uv)
{
    const uint32_t x = Min(static_cast<uint32_t>(uv.x * width), width - 1);
    const uint32_t y = Min(static_cast<uint32_t>(uv.y * height), height - 1);
    const float *row = conditionalCDF + size_t(y) * (size_t(width) + 1);
    const float rowProbability = marginalCDF[y + 1] - marginalCDF[y];
    const float columnProbability = row[x + 1] - row[x];
    return rowProbability * columnProbability * float(width) * float(height);
}

} // namespace ybi
