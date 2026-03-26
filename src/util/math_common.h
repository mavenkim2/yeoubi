#pragma once

#include "device/util.h"

namespace ybi
{

template <typename T>
YBI_DEVICE T Clamp(const T value, const T lo, const T hi)
{
    return value < lo ? lo : (value > hi ? hi : value);
}

YBI_DEVICE float Lerp(const float a, const float b, float t)
{
    return a * (1.0f - t) + b * t;
}

} // namespace ybi
