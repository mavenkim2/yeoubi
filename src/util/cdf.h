#pragma once

#include "util/array.h"
#include "util/assert.h"
#include "util/math_common.h"
#include "util/vec2.h"

#include <algorithm>
#include <cmath>

namespace ybi
{
namespace detail
{

inline void InitializeUniformCdf(int steps, float *cdf)
{
    YBI_ASSERT(cdf);
    YBI_ASSERT(steps > 0);

    const float invSteps = 1.0f / float(steps);
    for (int i = 0; i <= steps; ++i)
    {
        cdf[i] = float(i) * invSteps;
    }
}

inline void InitializeCdfFromWeights(const float *weights, int steps, float *cdf)
{
    YBI_ASSERT(weights);
    YBI_ASSERT(cdf);
    YBI_ASSERT(steps > 0);

    float totalWeight = 0.0f;
    for (int i = 0; i < steps; ++i)
    {
        totalWeight += std::max(weights[i], 0.0f);
    }

    if (totalWeight <= 0.0f)
    {
        InitializeUniformCdf(steps, cdf);
        return;
    }

    cdf[0] = 0.0f;
    float accumulated = 0.0f;
    for (int i = 0; i < steps; ++i)
    {
        accumulated += std::max(weights[i], 0.0f) / totalWeight;
        cdf[i + 1] = accumulated;
    }
    cdf[steps] = 1.0f;
}

} // namespace detail

template <typename Func>
void InitializeCDF1D(const Func &func, float start, float end, int steps, float *cdf)
{
    YBI_ASSERT(cdf);
    YBI_ASSERT(steps > 0);

    Array<float> weights{size_t(steps)};
    for (int i = 0; i < steps; ++i)
    {
        const float t = (float(i) + 0.5f) / float(steps);
        const float x = Lerp(start, end, t);
        weights[size_t(i)] = std::fabs(func(x));
    }

    detail::InitializeCdfFromWeights(weights.data(), steps, cdf);
}

inline void InitializeInverseCDF(
    const Array<float> &cdf, Array<float> &invCdf, float start, float end, int steps)
{
    YBI_ASSERT(steps > 0);
    YBI_ASSERT(cdf.size() >= size_t(steps + 1));

    invCdf.Resize(size_t(steps));
    const float invSteps = 1.0f / float(steps);
    const float *cdfBegin = cdf.begin();
    const float *cdfEnd = cdfBegin + steps + 1;

    for (int i = 0; i < steps; ++i)
    {
        const float value = (float(i) + 0.5f) * invSteps;
        const float *itr = std::upper_bound(cdfBegin, cdfEnd, value);
        int index = int(itr - cdfBegin) - 1;
        index = Clamp(index, 0, steps - 1);

        const float cdfLo = cdf[size_t(index)];
        const float cdfHi = cdf[size_t(index + 1)];
        const float cdfWidth = cdfHi - cdfLo;
        const float t = cdfWidth > 0.0f ? (value - cdfLo) / cdfWidth : 0.5f;
        invCdf[size_t(i)] = Lerp(start, end, (float(index) + t) * invSteps);
    }
}

template <typename Func>
void InitializeCDF2D(const Func &func,
                     Vec2 start,
                     Vec2 end,
                     int stepsU,
                     int stepsV,
                     float *conditional,
                     float *marginal)
{
    YBI_ASSERT(conditional);
    YBI_ASSERT(marginal);
    YBI_ASSERT(stepsU > 0);
    YBI_ASSERT(stepsV > 0);

    Array<float> rowWeights{size_t(stepsU)};
    Array<float> marginalWeights{size_t(stepsV)};

    for (int v = 0; v < stepsV; ++v)
    {
        const float tv = (float(v) + 0.5f) / float(stepsV);
        const float y = Lerp(start.y, end.y, tv);

        float rowIntegral = 0.0f;
        for (int u = 0; u < stepsU; ++u)
        {
            const float tu = (float(u) + 0.5f) / float(stepsU);
            const float x = Lerp(start.x, end.x, tu);
            const float weight = std::fabs(func(x, y));
            rowWeights[size_t(u)] = weight;
            rowIntegral += weight;
        }

        detail::InitializeCdfFromWeights(
            rowWeights.data(), stepsU, conditional + size_t(v) * size_t(stepsU + 1));
        marginalWeights[size_t(v)] = rowIntegral;
    }

    detail::InitializeCdfFromWeights(marginalWeights.data(), stepsV, marginal);
}

} // namespace ybi
