#pragma once

#include "util/math_common.h"

#include <cstddef>

namespace ybi
{
namespace render
{
namespace sampling
{

struct CdfSample1D
{
    int index = -1;
    float pmf = 0.0f;
    float remapped = 0.5f;
};

struct CdfSample2D
{
    int x = -1;
    int y = -1;
    float pmf = 0.0f;
    float u = 0.5f;
    float v = 0.5f;
};

YBI_DEVICE const float *CdfPtr(unsigned long long ptr)
{
    return ptr == 0ull ? nullptr : reinterpret_cast<const float *>(ptr);
}

YBI_DEVICE bool SampleCdf1D(const float *cdf, int count, float sample, CdfSample1D *outSample)
{
    if (!cdf || !outSample || count <= 0)
    {
        return false;
    }

    const float u = Clamp(sample, 0.0f, 0.99999994f);
    int begin = 0;
    int end = count;
    while (begin + 1 < end)
    {
        const int mid = begin + (end - begin) / 2;
        if (cdf[mid] <= u)
        {
            begin = mid;
        }
        else
        {
            end = mid;
        }
    }

    const int index = Clamp(begin, 0, count - 1);
    const float cdfLo = cdf[index];
    const float cdfHi = cdf[index + 1];
    const float pmf = cdfHi > cdfLo ? (cdfHi - cdfLo) : 0.0f;

    outSample->index = index;
    outSample->pmf = pmf;
    outSample->remapped = pmf > 0.0f ? Clamp((u - cdfLo) / pmf, 0.0f, 1.0f) : 0.5f;
    return true;
}

YBI_DEVICE float EvaluateCdf1DPmf(const float *cdf, int count, int index)
{
    if (!cdf || count <= 0)
    {
        return 0.0f;
    }
    const int i = Clamp(index, 0, count - 1);
    const float cdfLo = cdf[i];
    const float cdfHi = cdf[i + 1];
    return cdfHi > cdfLo ? (cdfHi - cdfLo) : 0.0f;
}

YBI_DEVICE bool SampleCdf2D(const float *conditional,
                            const float *marginal,
                            int width,
                            int height,
                            float sampleU,
                            float sampleV,
                            CdfSample2D *outSample)
{
    if (!conditional || !marginal || !outSample || width <= 0 || height <= 0)
    {
        return false;
    }

    CdfSample1D row = {};
    if (!SampleCdf1D(marginal, height, sampleV, &row))
    {
        return false;
    }

    const float *rowCdf = conditional + size_t(row.index) * size_t(width + 1);
    CdfSample1D column = {};
    if (!SampleCdf1D(rowCdf, width, sampleU, &column))
    {
        return false;
    }

    outSample->x = column.index;
    outSample->y = row.index;
    outSample->pmf = row.pmf * column.pmf;
    outSample->u = (float(column.index) + column.remapped) / float(width);
    outSample->v = (float(row.index) + row.remapped) / float(height);
    return true;
}

YBI_DEVICE float EvaluateCdf2DPmf(const float *conditional,
                                  const float *marginal,
                                  int width,
                                  int height,
                                  int x,
                                  int y)
{
    if (!conditional || !marginal || width <= 0 || height <= 0)
    {
        return 0.0f;
    }

    const int row = Clamp(y, 0, height - 1);
    const int col = Clamp(x, 0, width - 1);
    const float rowPmf = EvaluateCdf1DPmf(marginal, height, row);
    const float *rowCdf = conditional + size_t(row) * size_t(width + 1);
    const float columnPmf = EvaluateCdf1DPmf(rowCdf, width, col);
    return rowPmf * columnPmf;
}

} // namespace sampling
} // namespace render
} // namespace ybi
