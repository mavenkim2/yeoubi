#pragma once

#include "core/vec3f.h"
#include "device/cuda_compat.h"

namespace ybi
{

// Sources:
// Eric Heitz
// "Understanding the Masking-Shadowing Function in Microfacet-Based BRDFs", JCGT, 2014
YBI_CPU_GPU float EvaluateAnisoGGX(const Vec3f &v, const Vec2f &a)
{
    float x = (v.x / a.x) * (v.x / a.x);
    float y = (v.y / a.y) * (v.y / a.y);
    float z = v.z * v.z;
    float denom = x + y + z;
    denom *= denom;

    float result = 1.f / (YBI_PI * a.x * a.y * denom);
    return result;
}

YBI_CPU_GPU Vec3f SampleGGXVNDF(const Vec3f &wo, const Vec2f &a, const Vec2f &u)
{
    // fma(1 - u.y, 1 + wo.z, -wo.z)
    // z = 1 + wo.z - u.y - u.y * wo.z - wo.z = 1 - u.y - u.y * wo.z

    Vec3f woStd = Normalize(Vec3f(a * wo.Xy(), wo.z));

    float z = 1 - u.y - u.y * woStd.z;
    float sinTheta = sqrtf(Max(0.f, 1.f - z * z));

    float phi = 2 * YBI_PI * u.x;
    float cosPhi = Cos(phi);
    float sinPhi = Sin(phi);

    Vec3f result(sinTheta * cosPhi, sinTheta * sinPhi, z);
    result += woStd;

    result.x *= a.x;
    result.y *= a.y;
    return Normalize(result);
}

// See section 5.4 in Heitz, 2014
YBI_CPU_GPU float EvalAnisoSmithG1(const Vec3f &v, const Vec2f &alpha)
{
    const float zSqr = v.z * v.z;
    if (zSqr == 0.0f)
        return 0.0f;

    float x = alpha.x * v.x;
    x *= x;
    float y = alpha.y * v.y;
    y *= y;
    return 2.f / (1.f + sqrtf(1.0f + (x + y) / zSqr));
}

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

// NOTE: reflection only GGX
YBI_CPU_GPU bool SampleGGX(const Vec3f &viewDirection,
                           const Vec2f &a,
                           const Vec2f &u,
                           const float eta,
                           Vec3f &lightDirection,
                           Vec3f &weight)
{
    Vec3f halfway = SampleGGXVNDF(viewDirection, a, u);
    float hDotV = Dot(halfway, viewDirection);

    if (hDotV < 0.f)
    {
        return false;
    }

    lightDirection = Normalize(-viewDirection + 2.f * hDotV * halfway);

    // View/light direction must be in same hemisphere
    if (lightDirection.z * viewDirection.z <= 0.f)
    {
        return false;
    }

    float D = EvaluateAnisoGGX(halfway, a);

    if (D <= 0.f)
    {
        return false;
    }

    float G_masking = EvalAnisoSmithG1(viewDirection, a);
    float G_shadowing = EvalAnisoSmithG1(lightDirection, a);
    float F = FrDielectric(hDotV, eta);

    // weight = D * G * F / (4 * view.z * light.z) * light.z / pdf
    // pdf = D * G_masking * hdotv / (4.0f * hdotl * v.z);
    // hdotl = hdotv under reflection
    // pdf = D * G_masking / (4.f * v.z)

    // weight = G_shadowing * F
    // NOTE: uses uncorrelated Smith, which is G_masking * G_shadowing

    weight = Vec3f(G_shadowing * F);
    return true;
}

} // namespace ybi
