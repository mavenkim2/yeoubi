#pragma once

#include "device/cuda_compat.h"

namespace ybi
{

template <typename T>
YBI_CPU_GPU constexpr T Clamp(T value, T minValue, T maxValue)
{
    return value < minValue ? minValue : (value > maxValue ? maxValue : value);
}

template <typename T>
YBI_CPU_GPU constexpr T Min(T a, T b)
{
    return a < b ? a : b;
}

} // namespace ybi
