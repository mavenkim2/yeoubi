#include "render/integrator_light.h"

#include <cmath>
#include <cstdio>

namespace
{

using ybi::Vec3;
using ybi::render::integrator::DirectionToLatLongUv;

bool NearlyEqual(float actual, float expected, float eps)
{
    return std::fabs(actual - expected) <= eps;
}

bool NearlyEqualPeriodic(float actual, float expected, float period, float eps)
{
    const float delta = std::fabs(actual - expected);
    return delta <= eps || std::fabs(delta - period) <= eps;
}

bool ExpectUv(const char *label,
              const Vec3 &direction,
              float expectedU,
              bool periodicU,
              float expectedV,
              float eps)
{
    float u = 0.0f;
    float v = 0.0f;
    DirectionToLatLongUv(direction, &u, &v);

    const bool uOk = periodicU ? NearlyEqualPeriodic(u, expectedU, 1.0f, eps)
                               : NearlyEqual(u, expectedU, eps);
    const bool vOk = NearlyEqual(v, expectedV, eps);
    if (!uOk || !vOk)
    {
        std::fprintf(stderr,
                     "%s: got uv=(%.6f, %.6f) expected=(%.6f, %.6f)\n",
                     label,
                     u,
                     v,
                     expectedU,
                     expectedV);
        return false;
    }
    return true;
}

} // namespace

int main()
{
    const float eps = 1.0e-5f;
    bool ok = true;
    ok &= ExpectUv("+Z", Vec3(0.0f, 0.0f, 1.0f), 0.5f, false, 0.5f, eps);
    ok &= ExpectUv("+X", Vec3(1.0f, 0.0f, 0.0f), 0.25f, false, 0.5f, eps);
    ok &= ExpectUv("-Z", Vec3(0.0f, 0.0f, -1.0f), 0.0f, true, 0.5f, eps);
    ok &= ExpectUv("+Y", Vec3(0.0f, 1.0f, 0.0f), 0.5f, false, 0.0f, eps);
    ok &= ExpectUv("-Y", Vec3(0.0f, -1.0f, 0.0f), 0.5f, false, 1.0f, eps);
    if (!ok)
    {
        return 1;
    }

    std::puts("dome_light_uv_test: ok");
    return 0;
}
