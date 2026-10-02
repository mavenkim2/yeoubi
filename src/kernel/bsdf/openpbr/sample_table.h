#pragma once

#include "core/math.h"
#include "kernel/texture/sampling.h"

namespace ybi
{

constexpr int OpenPbrTableSize = 32;
constexpr int OpenPbrTableSizeMinusOne = OpenPbrTableSize - 1;
constexpr float OpenPbrIorMax = 2.5f;
constexpr float OpenPbrInverseIorMax = 1.f / OpenPbrIorMax;

YBI_GPU float OpenPbrIorToExactTableIndex(float ior)
{
    constexpr int halfTableSize = OpenPbrTableSize / 2;
    constexpr int halfTableSizeMinusOne = halfTableSize - 1;
    constexpr float inverseIorRange = 1.f / (OpenPbrIorMax - 1.f);

    if (ior < 1.f)
    {
        const float fraction = (1.f / ior - 1.f) * inverseIorRange;
        return halfTableSizeMinusOne - fraction * halfTableSizeMinusOne;
    }

    const float fraction = (ior - 1.f) * inverseIorRange;
    return halfTableSize + fraction * halfTableSizeMinusOne;
}

YBI_GPU float OpenPbrAlphaToExactTableIndex(float isotropicAlpha)
{
    return sqrtf(isotropicAlpha) * OpenPbrTableSizeMinusOne;
}

YBI_GPU float OpenPbrCosThetaToExactTableIndex(float cosTheta)
{
    return cosTheta * OpenPbrTableSizeMinusOne;
}

YBI_GPU float OpenPbrExactTableIndexToTextureCoordinate(float exactIndex)
{
    constexpr float inverseTableSize = 1.f / OpenPbrTableSize;
    constexpr float minCoordinate = 0.5f * inverseTableSize;
    constexpr float maxCoordinate = 1.f - minCoordinate;
    return Clamp(minCoordinate + exactIndex * inverseTableSize, minCoordinate, maxCoordinate);
}

YBI_GPU float ExtrapolateOpaqueDielectricEnergyComplement(float tableValue, float ior)
{
    if (ior > OpenPbrIorMax || ior < OpenPbrInverseIorMax)
    {
        constexpr float f0MaxRatio = (OpenPbrIorMax - 1.f) / (OpenPbrIorMax + 1.f);
        constexpr float f0Max = f0MaxRatio * f0MaxRatio;
        constexpr float inverseF0ExtrapolationRange = 1.f / (1.f - f0Max);
        const float f0 = Sqr((ior - 1.f) / (ior + 1.f));
        const float progressTowardZero = (f0 - f0Max) * inverseF0ExtrapolationRange;
        tableValue *= 1.f - progressTowardZero;
    }

    return tableValue;
}

YBI_GPU float SampleOpaqueDielectricTable(uint64_t texture,
                                          float ior,
                                          float cosTheta,
                                          float isotropicAlpha)
{
    const float iorCoordinate = OpenPbrExactTableIndexToTextureCoordinate(
        OpenPbrIorToExactTableIndex(ior));
    const float alphaCoordinate = OpenPbrExactTableIndexToTextureCoordinate(
        OpenPbrAlphaToExactTableIndex(isotropicAlpha));
    const float cosThetaCoordinate = OpenPbrExactTableIndexToTextureCoordinate(
        OpenPbrCosThetaToExactTableIndex(cosTheta));

    const float tableValue = SampleDeviceTexture3D(texture,
                                                   TextureFormat::R16_UNORM,
                                                   cosThetaCoordinate,
                                                   alphaCoordinate,
                                                   iorCoordinate)
                                 .x;
    return ExtrapolateOpaqueDielectricEnergyComplement(tableValue, ior);
}

YBI_GPU float SampleOpaqueDielectricAverageTable(uint64_t texture,
                                                 float ior,
                                                 float isotropicAlpha)
{
    const float iorCoordinate = OpenPbrExactTableIndexToTextureCoordinate(
        OpenPbrIorToExactTableIndex(ior));
    const float alphaCoordinate = OpenPbrExactTableIndexToTextureCoordinate(
        OpenPbrAlphaToExactTableIndex(isotropicAlpha));

    const float tableValue =
        SampleDeviceTexture(texture, TextureFormat::R16_UNORM, alphaCoordinate, iorCoordinate, 0).x;
    return ExtrapolateOpaqueDielectricEnergyComplement(tableValue, ior);
}

} // namespace ybi
