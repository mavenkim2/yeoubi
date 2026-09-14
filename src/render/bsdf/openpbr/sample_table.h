#pragma once

#include "core/math.h"
#include "render/texture_sampling.h"

namespace ybi
{

constexpr int OpenPbrTableSize = 32;
constexpr int OpenPbrTableSizeMinusOne = OpenPbrTableSize - 1;
constexpr float OpenPbrIorMax = 2.5f;
constexpr float OpenPbrInverseIorMax = 1.f / OpenPbrIorMax;

YBI_GPU float SampleOpaqueDielectricTable(uint64_t texture,
                                          float ior,
                                          float cosTheta,
                                          float isotropicAlpha)
{
    constexpr int halfTableSize = OpenPbrTableSize / 2;
    constexpr int halfTableSizeMinusOne = halfTableSize - 1;
    constexpr float inverseIorRange = 1.f / (OpenPbrIorMax - 1.f);

    float exactIorIndex;
    if (ior < 1.f)
    {
        float fraction = (1.f / ior - 1.f) * inverseIorRange;
        exactIorIndex = halfTableSizeMinusOne - fraction * halfTableSizeMinusOne;
    }
    else
    {
        float fraction = (ior - 1.f) * inverseIorRange;
        exactIorIndex = halfTableSize + fraction * halfTableSizeMinusOne;
    }

    float exactAlphaIndex = sqrtf(isotropicAlpha) * OpenPbrTableSizeMinusOne;
    float exactCosThetaIndex = cosTheta * OpenPbrTableSizeMinusOne;

    constexpr float inverseTableSize = 1.f / OpenPbrTableSize;
    constexpr float minCoordinate = 0.5f * inverseTableSize;
    constexpr float maxCoordinate = 1.f - minCoordinate;
    float iorCoordinate =
        Clamp(minCoordinate + exactIorIndex * inverseTableSize, minCoordinate, maxCoordinate);
    float alphaCoordinate =
        Clamp(minCoordinate + exactAlphaIndex * inverseTableSize, minCoordinate, maxCoordinate);
    float cosThetaCoordinate = Clamp(
        minCoordinate + exactCosThetaIndex * inverseTableSize, minCoordinate, maxCoordinate);

    float tableValue = SampleDeviceTexture3D(texture,
                                             TextureFormat::R16_UNORM,
                                             cosThetaCoordinate,
                                             alphaCoordinate,
                                             iorCoordinate)
                           .x;

    if (ior > OpenPbrIorMax || ior < OpenPbrInverseIorMax)
    {
        constexpr float f0MaxRatio = (OpenPbrIorMax - 1.f) / (OpenPbrIorMax + 1.f);
        constexpr float f0Max = f0MaxRatio * f0MaxRatio;
        constexpr float inverseF0ExtrapolationRange = 1.f / (1.f - f0Max);
        float f0 = Sqr((ior - 1.f) / (ior + 1.f));
        float progressTowardZero = (f0 - f0Max) * inverseF0ExtrapolationRange;
        tableValue *= 1.f - progressTowardZero;
    }

    return tableValue;
}

} // namespace ybi
