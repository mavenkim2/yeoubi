#pragma once

#include "core/vec3f.h"

namespace ybi 
{

YBI_CPU_GPU float FrDielectric(float cosTheta_i, float eta)
{
    cosTheta_i = Clamp(cosTheta_i, -1.f, 1.f);

    if (cosTheta_i < 0)
    {
        cosTheta_i = -cosTheta_i;
        eta = 1.f / eta;
    }
    float sin2Theta_i = Max(0.f, 1 - cosTheta_i * cosTheta_i);
    float sin2Theta_t = sin2Theta_i / (eta * eta);

    // TIR
    if (sin2Theta_t >= 1.f)
    {
        return 1.f;
    }

    float cosTheta_t = sqrtf(1 - sin2Theta_t);
    float rParallel = (eta * cosTheta_i - cosTheta_t) / (eta * cosTheta_i + cosTheta_t);
    float rPerp = (cosTheta_i - eta * cosTheta_t) / (cosTheta_i + eta * cosTheta_t);
    float Fr = 0.5f * (rPerp * rPerp + rParallel * rParallel);
    return Fr;
}

YBI_GPU Vec2f ComputeAnisotropicRoughness(float roughness, float anisotropy)
{
    Vec2f anisotropicRoughness;
    anisotropicRoughness.x = Sqr(roughness) * sqrtf(2.f / (1.f + Sqr(1.0f - anisotropy)));
    anisotropicRoughness.y = (1.0f - anisotropy) * anisotropicRoughness.x;

    // TODO: support dirac delta lobes
    anisotropicRoughness = Max(anisotropicRoughness, Sqr(0.001f));
    return anisotropicRoughness;
}

YBI_GPU float CalculateThinWallReflectance(float cosTheta, float reflectionIor)
{
    const float frontReflectance = FrDielectric(Abs(cosTheta), reflectionIor);
    return 2.f * frontReflectance / (1.f + frontReflectance);
}

YBI_GPU Vec3f CalculateThinWallTransmissionColor(const Vec3f &viewDirection,
                                                 float specularIor,
                                                 const Vec3f &transmissionColor)
{
    float cosTheta = Abs(viewDirection.z);
    float sin2Theta_i = 1.f - Sqr(cosTheta);

    // TODO: this ior has to change when the coat layer is added
    float sin2Theta_t = sin2Theta_i / Sqr(specularIor);

    if (sin2Theta_t >= 1.f)
    {
        return Vec3f(0.f);
    }

    float cos_t = sqrtf(1.f - sin2Theta_t);
    float pathLength = 1.f / cos_t;
    Vec3f thinWallTransmissionColor = Pow(transmissionColor, pathLength);
    return thinWallTransmissionColor;
}

// https://blog.selfshadow.com/publications/s2017-shading-course/imageworks/s2017_pbs_imageworks_slides_v2.pdf,
// slide 18.
YBI_GPU float ApproximateAverageFresnel(const float eta)
{
    if (eta > 1.0f)
        return (eta - 1.0f) / (4.08567f + 1.00071f * eta);
    else
        return 0.997118f + 0.1014f * eta - 0.965241f * (eta * eta) - 0.130607f * (eta * eta) * eta;
}

// NOTE: ior weighting is only used for reflection, it is not applied to refraction
YBI_GPU float ApplyWeightToIor(float ior, float weight)
{
    float fNormalIncidence = Sqr((1 - ior) / (1 + ior));
    float eps = sqrtf(Min(weight * fNormalIncidence, 0.9999f));
    eps = ior >= 1.f ? eps : -eps;
    float result = (1 + eps) / (1 - eps);
    return result;
}

}
