#pragma once

#include "device/cuda_compat.h"
#include <cmath>

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

template <typename T>
YBI_CPU_GPU constexpr T Max(T a, T b)
{
    return a > b ? a : b;
}

YBI_CPU_GPU float Floor(float a)
{
    return floorf(a);
}

YBI_CPU_GPU float Abs(float a)
{
    return a < 0.f ? -a : a;
}

YBI_CPU_GPU constexpr float DegreesToRadians(float degrees)
{
    return degrees * 0.017453292519943295769f;
}

template <typename T>
YBI_CPU_GPU void Swap(T &a, T &b)
{
    T temp = a;
    a = b;
    b = temp;
}

} // namespace ybi
