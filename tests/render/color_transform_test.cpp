#include "render/color_transform.h"

#include <cmath>
#include <cstdio>

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
    ok &= ExpectFloat("LinearToSrgb(0.0)", ybi::render::LinearToSrgb(0.0f), 0.0f, 1.0e-7f);
    ok &= ExpectFloat("LinearToSrgb(0.0031308)",
                      ybi::render::LinearToSrgb(0.0031308f),
                      0.04044994f,
                      1.0e-6f);
    ok &= ExpectFloat(
        "LinearToSrgb(0.18)", ybi::render::LinearToSrgb(0.18f), 0.46135613f, 1.0e-6f);

    const ybi::Vec3 zero = ybi::render::DisplayMapPathRadiance(ybi::Vec3(0.0f, 0.0f, 0.0f));
    ok &= ExpectFloat("DisplayMapPathRadiance(0).x", zero.x, 0.0f, 1.0e-7f);
    ok &= ExpectFloat("DisplayMapPathRadiance(0).y", zero.y, 0.0f, 1.0e-7f);
    ok &= ExpectFloat("DisplayMapPathRadiance(0).z", zero.z, 0.0f, 1.0e-7f);

    const float dark = 0.0f;
    const float bright = 8.0f;
    const float mappedAverage =
        ybi::render::DisplayMapPathRadiance(ybi::Vec3((dark + bright) * 0.5f, 0.0f, 0.0f)).x;
    const float averagedMapped =
        0.5f *
        (ybi::render::DisplayMapPathRadiance(ybi::Vec3(dark, 0.0f, 0.0f)).x +
         ybi::render::DisplayMapPathRadiance(ybi::Vec3(bright, 0.0f, 0.0f)).x);
    ok &= ExpectGreater("Average-before-display", mappedAverage, averagedMapped, 0.25f);

    if (!ok)
    {
        return 1;
    }

    std::puts("color_transform_test: ok");
    return 0;
}
