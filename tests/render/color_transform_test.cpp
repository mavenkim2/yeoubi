#include "render/color_transform.h"

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>

namespace
{

bool NearlyEqual(float actual, float expected, float eps)
{
    return std::fabs(actual - expected) <= eps;
}

bool ExpectFloat(const char *label, float actual, float expected, float eps)
{
    if (NearlyEqual(actual, expected, eps))
    {
        return true;
    }

    std::fprintf(stderr, "%s: got %.7f expected %.7f\n", label, actual, expected);
    return false;
}

bool ExpectGreater(const char *label, float lhs, float rhs, float margin)
{
    if (lhs > rhs + margin)
    {
        return true;
    }

    std::fprintf(stderr,
                 "%s: expected %.7f > %.7f by at least %.7f\n",
                 label,
                 lhs,
                 rhs,
                 margin);
    return false;
}

} // namespace

int main()
{
    bool ok = true;
    const std::vector<float> linearRgb = {
        0.0f, 0.0f, 0.0f,
        0.18f, 0.18f, 0.18f,
        4.0f, 4.0f, 4.0f,
    };
    std::vector<uint8_t> rgba;
    std::string error;
    ok &= ybi::render::ApplyAcesSdrDisplayTransform(linearRgb, 3, 1, &rgba, &error);
    if (!ok)
    {
        std::fprintf(stderr, "ApplyAcesSdrDisplayTransform failed: %s\n", error.c_str());
    }
    else
    {
        ok &= ExpectFloat("black.r", static_cast<float>(rgba[0]), 0.0f, 0.5f);
        ok &= ExpectFloat("black.g", static_cast<float>(rgba[1]), 0.0f, 0.5f);
        ok &= ExpectFloat("black.b", static_cast<float>(rgba[2]), 0.0f, 0.5f);

        ok &= ExpectFloat("middle neutral rg",
                          static_cast<float>(rgba[4]),
                          static_cast<float>(rgba[5]),
                          0.5f);
        ok &= ExpectFloat("middle neutral gb",
                          static_cast<float>(rgba[5]),
                          static_cast<float>(rgba[6]),
                          0.5f);
        ok &= ExpectFloat("bright neutral rg",
                          static_cast<float>(rgba[8]),
                          static_cast<float>(rgba[9]),
                          0.5f);
        ok &= ExpectFloat("bright neutral gb",
                          static_cast<float>(rgba[9]),
                          static_cast<float>(rgba[10]),
                          0.5f);

        ok &= ExpectGreater("middle > black",
                            static_cast<float>(rgba[4]),
                            static_cast<float>(rgba[0]),
                            1.0f);
        ok &= ExpectGreater("bright > middle",
                            static_cast<float>(rgba[8]),
                            static_cast<float>(rgba[4]),
                            1.0f);
    }

    if (!ok)
    {
        return 1;
    }

    std::puts("color_transform_test: ok");
    return 0;
}
