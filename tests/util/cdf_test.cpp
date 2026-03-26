#include "util/array.h"
#include "util/cdf.h"

#include <cmath>
#include <cstdio>

namespace
{

bool NearlyEqual(float actual, float expected, float eps = 1.0e-6f)
{
    return std::fabs(actual - expected) <= eps;
}

bool TestInitializeCDF1DConstant()
{
    float cdf[5] = {};
    ybi::InitializeCDF1D([](float) { return 1.0f; }, 0.0f, 1.0f, 4, cdf);

    const float expected[5] = {0.0f, 0.25f, 0.5f, 0.75f, 1.0f};
    for (int i = 0; i < 5; ++i)
    {
        if (!NearlyEqual(cdf[i], expected[i]))
        {
            std::fprintf(stderr, "cdf_test: constant 1D mismatch at %d: %.6f\n", i, cdf[i]);
            return false;
        }
    }
    return true;
}

bool TestInitializeCDF1DZeroFallback()
{
    float cdf[4] = {};
    ybi::InitializeCDF1D([](float) { return 0.0f; }, -1.0f, 1.0f, 3, cdf);

    const float expected[4] = {0.0f, 1.0f / 3.0f, 2.0f / 3.0f, 1.0f};
    for (int i = 0; i < 4; ++i)
    {
        if (!NearlyEqual(cdf[i], expected[i]))
        {
            std::fprintf(stderr, "cdf_test: zero fallback mismatch at %d: %.6f\n", i, cdf[i]);
            return false;
        }
    }
    return true;
}

bool TestInitializeInverseCDFUniform()
{
    ybi::Array<float> cdf(5);
    cdf[0] = 0.0f;
    cdf[1] = 0.25f;
    cdf[2] = 0.5f;
    cdf[3] = 0.75f;
    cdf[4] = 1.0f;

    ybi::Array<float> inv;
    ybi::InitializeInverseCDF(cdf, inv, 0.0f, 1.0f, 4);

    const float expected[4] = {0.125f, 0.375f, 0.625f, 0.875f};
    for (int i = 0; i < 4; ++i)
    {
        if (!NearlyEqual(inv[size_t(i)], expected[i]))
        {
            std::fprintf(stderr, "cdf_test: inverse mismatch at %d: %.6f\n", i, inv[size_t(i)]);
            return false;
        }
    }
    return true;
}

bool TestInitializeCDF2DConstant()
{
    float conditional[2 * 3] = {};
    float marginal[3] = {};
    ybi::InitializeCDF2D(
        [](float, float) { return 1.0f; }, ybi::Vec2(0.0f), ybi::Vec2(1.0f), 2, 2, conditional, marginal);

    const float expectedRow[3] = {0.0f, 0.5f, 1.0f};
    for (int row = 0; row < 2; ++row)
    {
        for (int i = 0; i < 3; ++i)
        {
            const float actual = conditional[row * 3 + i];
            if (!NearlyEqual(actual, expectedRow[i]))
            {
                std::fprintf(stderr,
                             "cdf_test: 2D conditional mismatch row=%d idx=%d val=%.6f\n",
                             row,
                             i,
                             actual);
                return false;
            }
        }
    }

    const float expectedMarginal[3] = {0.0f, 0.5f, 1.0f};
    for (int i = 0; i < 3; ++i)
    {
        if (!NearlyEqual(marginal[i], expectedMarginal[i]))
        {
            std::fprintf(stderr, "cdf_test: 2D marginal mismatch at %d: %.6f\n", i, marginal[i]);
            return false;
        }
    }
    return true;
}

bool TestInitializeCDF2DZeroRowFallback()
{
    float conditional[2 * 3] = {};
    float marginal[3] = {};
    ybi::InitializeCDF2D([](float, float y) { return y > 0.5f ? 1.0f : 0.0f; },
                         ybi::Vec2(0.0f),
                         ybi::Vec2(1.0f),
                         2,
                         2,
                         conditional,
                         marginal);

    const float expectedConditional[6] = {0.0f, 0.5f, 1.0f, 0.0f, 0.5f, 1.0f};
    for (int i = 0; i < 6; ++i)
    {
        if (!NearlyEqual(conditional[i], expectedConditional[i]))
        {
            std::fprintf(stderr, "cdf_test: zero-row conditional mismatch at %d: %.6f\n", i, conditional[i]);
            return false;
        }
    }

    const float expectedMarginal[3] = {0.0f, 0.0f, 1.0f};
    for (int i = 0; i < 3; ++i)
    {
        if (!NearlyEqual(marginal[i], expectedMarginal[i]))
        {
            std::fprintf(stderr, "cdf_test: zero-row marginal mismatch at %d: %.6f\n", i, marginal[i]);
            return false;
        }
    }
    return true;
}

} // namespace

int main()
{
    bool ok = true;
    ok &= TestInitializeCDF1DConstant();
    ok &= TestInitializeCDF1DZeroFallback();
    ok &= TestInitializeInverseCDFUniform();
    ok &= TestInitializeCDF2DConstant();
    ok &= TestInitializeCDF2DZeroRowFallback();
    if (!ok)
    {
        return 1;
    }

    std::puts("cdf_test: ok");
    return 0;
}
