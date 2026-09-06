// EON implementation adapted from https://github.com/portsmouth/EON-diffuse.
// Copyright (c) 2024 Jamie Portsmouth
// SPDX-License-Identifier: MIT
// Full license: licenses/EON-diffuse.txt (relative to the repository root).
#pragma once

#include "core/math.h"
#include "render/sampling.h"
#include <cmath>
namespace ybi
{

// Jamie Portsmouth, Peter Kutz, and Stephen Hill.
// "EON: A Practical Energy-Preserving Rough Diffuse BRDF", JCGT, 2025.
// https://jcgt.org/published/0014/01/06/
static inline constexpr float OREN_NAYAR_CONSTANT_1 = 0.5f - 2.f / (3.f * YBI_PI);
static inline constexpr float OREN_NAYAR_CONSTANT_2 = 2.f / 3.f - 28.f / (15.f * YBI_PI);

YBI_CPU_GPU float EON_exact(float cosTheta, float r)
{
    float Af = 1.f / (1.f + OREN_NAYAR_CONSTANT_1 * r);
    float Bf = r * Af;

    cosTheta = Clamp(cosTheta, 0.f, 1.f);
    float sinTheta = sqrtf(1.f - cosTheta * cosTheta);

    // Uses cancellation free version of 1 - sin^3
    // 1 - sin^3 = (1 - s)(1 + s + s^2)
    // mu^2 = 1 - s^2 = (1 + s)(1 - s) -> (1-s) = mu^2 / (1 + s)
    float G = sinTheta * (acosf(cosTheta) - sinTheta * cosTheta) +
              (2.f / 3.f) * (sinTheta * cosTheta * (1.f + sinTheta + sinTheta * sinTheta) /
                                 (1.f + sinTheta) -
                             sinTheta);

    float result = Af + Bf * YBI_RCP_PI * G;
    return result;
}

YBI_CPU_GPU Vec3f EvaluateOrenNayar(const Vec3f &viewDirection,
                                    const Vec3f &lightDirection,
                                    const Vec3f &rho,
                                    float r)
{
    float cosThetaV = viewDirection.z;
    float cosThetaL = lightDirection.z;

    float s = Dot(viewDirection, lightDirection) - cosThetaV * cosThetaL;
    float sOverF = s > 0.f ? s / Max(cosThetaV, cosThetaL) : s;

    float Af = 1.f / (1.f + OREN_NAYAR_CONSTANT_1 * r);
    float Bf = r * Af;

    Vec3f f_ss = rho * YBI_RCP_PI * (Af + Bf * sOverF);
    float EFo = EON_exact(cosThetaV, r);
    float EFi = EON_exact(cosThetaL, r);
    float avgEF = Af + Bf * OREN_NAYAR_CONSTANT_2;

    Vec3f rho_ms = (rho * rho) * avgEF / (1.f - rho * (1.f - avgEF));

    constexpr float eps = 1.0e-7f;

    float oneMinusEFo = Max(1 - EFo, eps);
    float oneMinusEFi = Max(1 - EFi, eps);
    float oneMinusAvgEF = Max(1 - avgEF, eps);

    Vec3f f_ms = rho_ms * YBI_RCP_PI * oneMinusEFo * oneMinusEFi / oneMinusAvgEF;
    Vec3f result = f_ss + f_ms;
    return result;
}

YBI_CPU_GPU bool SampleOrenNayar(const Vec3f &viewDirection,
                                 const Vec3f &baseColor,
                                 const Vec2f &u,
                                 Vec3f &lightDirection,
                                 Vec3f &weight)
{
    Vec3f lightDirectionLocal = SampleCosineHemisphere(u);
    weight = EvaluateOrenNayar(viewDirection, lightDirectionLocal, baseColor, 0.0f) * YBI_PI;

    lightDirection = lightDirectionLocal;
    return true;
}

} // namespace ybi
