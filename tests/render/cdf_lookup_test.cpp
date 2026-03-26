#include "render/sampling/cdf_lookup.h"
#include "util/cdf.h"

#include <cmath>
#include <cstdio>

namespace
{

using ybi::render::sampling::CdfSample1D;
using ybi::render::sampling::CdfSample2D;
using ybi::render::sampling::EvaluateCdf2DPmf;
using ybi::render::sampling::SampleCdf1D;
using ybi::render::sampling::SampleCdf2D;

bool NearlyEqual(float actual, float expected, float eps = 1.0e-6f)
{
    return std::fabs(actual - expected) <= eps;
}

bool TestSampleCdf1D()
{
    float cdf[4] = {0.0f, 0.25f, 0.75f, 1.0f};
    CdfSample1D sample = {};
    if (!SampleCdf1D(cdf, 3, 0.50f, &sample))
    {
        std::fprintf(stderr, "cdf_lookup_test: SampleCdf1D failed\n");
        return false;
    }
    if (sample.index != 1 || !NearlyEqual(sample.pmf, 0.50f) ||
        !NearlyEqual(sample.remapped, 0.50f))
    {
        std::fprintf(stderr,
                     "cdf_lookup_test: bad 1D sample index=%d pmf=%.6f remapped=%.6f\n",
                     sample.index,
                     sample.pmf,
                     sample.remapped);
        return false;
    }
    return true;
}

bool TestSampleCdf2D()
{
    float conditional[2 * 3] = {};
    float marginal[3] = {};
    ybi::InitializeCDF2D([](float u, float v) { return (u < 0.5f ? 1.0f : 3.0f) * (v < 0.5f ? 2.0f : 1.0f); },
                         ybi::Vec2(0.0f),
                         ybi::Vec2(1.0f),
                         2,
                         2,
                         conditional,
                         marginal);

    CdfSample2D sample = {};
    if (!SampleCdf2D(conditional, marginal, 2, 2, 0.80f, 0.20f, &sample))
    {
        std::fprintf(stderr, "cdf_lookup_test: SampleCdf2D failed\n");
        return false;
    }
    if (sample.x != 1 || sample.y != 0)
    {
        std::fprintf(stderr,
                     "cdf_lookup_test: bad 2D sample xy=(%d,%d)\n",
                     sample.x,
                     sample.y);
        return false;
    }

    const float pmf00 = EvaluateCdf2DPmf(conditional, marginal, 2, 2, 0, 0);
    const float pmf10 = EvaluateCdf2DPmf(conditional, marginal, 2, 2, 1, 0);
    const float pmf01 = EvaluateCdf2DPmf(conditional, marginal, 2, 2, 0, 1);
    const float pmf11 = EvaluateCdf2DPmf(conditional, marginal, 2, 2, 1, 1);
    if (!NearlyEqual(pmf00 + pmf10 + pmf01 + pmf11, 1.0f) ||
        !(pmf10 > pmf00 && pmf00 > pmf01 && pmf11 > pmf01))
    {
        std::fprintf(stderr,
                     "cdf_lookup_test: bad PMFs %.6f %.6f %.6f %.6f\n",
                     pmf00,
                     pmf10,
                     pmf01,
                     pmf11);
        return false;
    }
    return true;
}

} // namespace

int main()
{
    bool ok = true;
    ok &= TestSampleCdf1D();
    ok &= TestSampleCdf2D();
    if (!ok)
    {
        return 1;
    }

    std::puts("cdf_lookup_test: ok");
    return 0;
}
