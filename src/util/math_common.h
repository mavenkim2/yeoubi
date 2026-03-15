#pragma once

#include "util/forceinline.h"

#ifndef YBI_DEVICE
#if defined(__CUDACC__)
#define YBI_DEVICE __host__ __device__ __forceinline__
#define YBI_DEVICE_FORCEINLINE __device__ __forceinline__
#else
#define YBI_DEVICE __forceinline
#define YBI_DEVICE_FORCEINLINE __forceinline
#endif
#endif

namespace ybi
{

YBI_DEVICE float Lerp(const float a, const float b, float t)
{
    return a * (1.0f - t) + b * t;
}

} // namespace ybi
