#include "render/sampling/cdf.h"

#include "core/assert.h"
#include <atomic>
#include <cmath>
#include <tbb/parallel_for.h>

namespace ybi
{

float CreatePiecewise1DCDF(const float *values, uint32_t numValues, float *cdf)
{
    YBI_ASSERT(numValues > 0);

    float total = 0.f;
    for (uint32_t valueIndex = 0; valueIndex < numValues; valueIndex++)
        total += values[valueIndex];

    float cumulative = 0.f;
    for (uint32_t valueIndex = 0; valueIndex < numValues; valueIndex++)
    {
        cdf[valueIndex] = total > 0.f ? cumulative / total
                                    : float(valueIndex) / float(numValues);
        cumulative += values[valueIndex];
    }
    cdf[numValues] = 1.f;
    return total;
}

bool CreatePiecewise2DCDF(const std::vector<float> &values,
                          uint32_t width,
                          uint32_t height,
                          std::vector<float> &conditionalCDF,
                          std::vector<float> &marginalCDF,
                          float &totalWeight)
{
    if (width == 0 || height == 0 || values.size() != size_t(width) * height)
        return false;

    const size_t cdfWidth = size_t(width) + 1;
    conditionalCDF.resize(cdfWidth * height);
    marginalCDF.resize(size_t(height) + 1);
    std::vector<float> rowTotals(height);
    std::atomic<bool> valid{true};

    tbb::parallel_for(size_t(0), size_t(height), [&](size_t y) {
        const float *row = values.data() + y * width;
        for (uint32_t x = 0; x < width; x++)
        {
            if (!std::isfinite(row[x]) || row[x] < 0.f)
            {
                valid.store(false, std::memory_order_relaxed);
                return;
            }
        }
        rowTotals[y] = CreatePiecewise1DCDF(row, width, conditionalCDF.data() + y * cdfWidth);
        if (!std::isfinite(rowTotals[y]))
            valid.store(false, std::memory_order_relaxed);
    });

    if (!valid.load(std::memory_order_relaxed))
        return false;

    totalWeight = CreatePiecewise1DCDF(rowTotals.data(), height, marginalCDF.data());
    return std::isfinite(totalWeight);
}

} // namespace ybi
