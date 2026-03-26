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
