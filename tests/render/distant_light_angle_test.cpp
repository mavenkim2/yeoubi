#include "render/integrator_light.h"

#include <cmath>
#include <cstdio>

namespace
{

using ybi::Float3x4;
using ybi::LaunchParams;
using ybi::PackedLight;
using ybi::Vec3;
using ybi::render::integrator::DirectLightSample;
using ybi::render::integrator::DirectionInsideDistantLightCap;
using ybi::render::integrator::DistantLightDirectionalPdf;
using ybi::render::integrator::EvaluateDistantLightRadiance;
using ybi::render::integrator::SampleDirectLight;

bool NearlyEqual(float actual, float expected, float eps)
{
    return std::fabs(actual - expected) <= eps;
}

PackedLight MakeDistantLight(float angleDegrees)
{
    PackedLight light = {};
    light.type = static_cast<unsigned int>(ybi::LightType::Distant);
    light.worldFromLocal = Float3x4::Identity();
    light.emissionScale = 3.0f;
    light.color = Vec3(2.0f, 1.0f, 0.5f);
    light.selectionWeight = 1.0f;

    const float angleRadians = angleDegrees * (ybi::kPi / 180.0f);
    light.angleRadians = angleRadians;
    if (angleRadians <= 1.0e-6f)
    {
        light.cosThetaMax = 1.0f;
        light.solidAngle = 0.0f;
        return light;
    }

    const float thetaMax = 0.5f * angleRadians;
    light.cosThetaMax = std::cos(thetaMax);
    light.solidAngle = ybi::kTwoPi * (1.0f - light.cosThetaMax);
    return light;
}

Vec3 DirectionAtAngleDegrees(float degrees)
{
    const float radians = degrees * (ybi::kPi / 180.0f);
    return Vec3(std::sin(radians), 0.0f, std::cos(radians));
}

bool TestCapMembership()
{
    const PackedLight light = MakeDistantLight(20.0f);
    const Vec3 inside = DirectionAtAngleDegrees(9.5f);
    const Vec3 outside = DirectionAtAngleDegrees(10.5f);
    const Vec3 axis(0.0f, 0.0f, 1.0f);

    if (!DirectionInsideDistantLightCap(light, axis) ||
        !DirectionInsideDistantLightCap(light, inside) ||
        DirectionInsideDistantLightCap(light, outside))
    {
        std::fprintf(stderr, "distant_light_angle_test: cap membership failed\n");
        return false;
    }
    return true;
}

bool TestFiniteAngleSamplePdf()
{
    PackedLight light = MakeDistantLight(20.0f);
    LaunchParams params = {};
    params.lights = reinterpret_cast<unsigned long long>(&light);
    params.lightCount = 1;

    const float expectedPdf = 1.0f / light.solidAngle;
    const float helperPdf = DistantLightDirectionalPdf(light);
    if (!NearlyEqual(helperPdf, expectedPdf, 1.0e-6f))
    {
        std::fprintf(stderr,
                     "distant_light_angle_test: helper pdf %.8f expected %.8f\n",
                     helperPdf,
                     expectedPdf);
        return false;
    }

    unsigned int rngState = 1u;
    for (int i = 0; i < 64; ++i)
    {
        DirectLightSample sample = {};
        if (!SampleDirectLight(params, Vec3(0.0f, 0.0f, 0.0f), rngState, &sample))
        {
            std::fprintf(stderr, "distant_light_angle_test: failed finite-angle sample\n");
            return false;
        }
        if (sample.isDeltaLight || !DirectionInsideDistantLightCap(light, sample.wi) ||
            !NearlyEqual(sample.pdf, expectedPdf, 1.0e-5f))
        {
            std::fprintf(stderr,
                         "distant_light_angle_test: bad finite sample delta=%d pdf=%.8f\n",
                         sample.isDeltaLight ? 1 : 0,
                         sample.pdf);
            return false;
        }
    }

    const Vec3 insideRadiance = EvaluateDistantLightRadiance(light, Vec3(0.0f, 0.0f, 1.0f));
    const Vec3 outsideRadiance = EvaluateDistantLightRadiance(light, DirectionAtAngleDegrees(25.0f));
    if (!NearlyEqual(insideRadiance.x, 6.0f, 1.0e-6f) ||
        !NearlyEqual(outsideRadiance.x, 0.0f, 1.0e-6f))
    {
        std::fprintf(stderr, "distant_light_angle_test: bad miss radiance evaluation\n");
        return false;
    }
    return true;
}

bool TestZeroAnglePreservesDelta()
{
    PackedLight light = MakeDistantLight(0.0f);
    LaunchParams params = {};
    params.lights = reinterpret_cast<unsigned long long>(&light);
    params.lightCount = 1;

    unsigned int rngState = 7u;
    DirectLightSample sample = {};
    if (!SampleDirectLight(params, Vec3(0.0f, 0.0f, 0.0f), rngState, &sample))
    {
        std::fprintf(stderr, "distant_light_angle_test: failed zero-angle sample\n");
        return false;
    }
    if (!sample.isDeltaLight || !NearlyEqual(sample.pdf, 1.0f, 1.0e-6f) ||
        !NearlyEqual(sample.wi.z, 1.0f, 1.0e-6f))
    {
        std::fprintf(stderr,
                     "distant_light_angle_test: zero-angle delta=%d pdf=%.8f wi=(%.6f %.6f %.6f)\n",
                     sample.isDeltaLight ? 1 : 0,
                     sample.pdf,
                     sample.wi.x,
                     sample.wi.y,
                     sample.wi.z);
        return false;
    }
    return true;
}

} // namespace

int main()
{
    bool ok = true;
    ok &= TestCapMembership();
    ok &= TestFiniteAngleSamplePdf();
    ok &= TestZeroAnglePreservesDelta();
    if (!ok)
    {
        return 1;
    }

    std::puts("distant_light_angle_test: ok");
    return 0;
}
