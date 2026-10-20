#include "render/sampling/cdf.h"

#include "core/assert.h"
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

float CreatePiecewise2DCDF(const std::vector<float> &values,
                          uint32_t width,
                          uint32_t height,
                          std::vector<float> &conditionalCDF,
                          std::vector<float> &marginalCDF)
{
    YBI_ASSERT(width > 0 && height > 0);
    YBI_ASSERT(values.size() == size_t(width) * height);

    const size_t cdfWidth = size_t(width) + 1;
    conditionalCDF.resize(cdfWidth * height);
    marginalCDF.resize(size_t(height) + 1);
    std::vector<float> rowTotals(height);

    tbb::parallel_for(size_t(0), size_t(height), [&](size_t y) {
        const float *row = values.data() + y * width;
        rowTotals[y] = CreatePiecewise1DCDF(row, width, conditionalCDF.data() + y * cdfWidth);
    });

    return CreatePiecewise1DCDF(rowTotals.data(), height, marginalCDF.data());
}

} // namespace ybi
