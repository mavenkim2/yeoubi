#pragma once

#include "util/forceinline.h"

#ifndef YBI_INTEGRATOR_HD
#if defined(__CUDACC__)
#define YBI_INTEGRATOR_HD __host__ __device__ __forceinline__
#else
#define YBI_INTEGRATOR_HD __forceinline
#endif
#endif
