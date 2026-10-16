#pragma once

#include "device/cuda_compat.h"

#if defined(_MSC_VER) && !defined(__CUDACC_RTC__)
#include <intrin.h>
#endif

namespace ybi
{

YBI_CPU_GPU uint32_t AtomicOr(uint32_t *destination, uint32_t value)
{
#if defined(__CUDA_ARCH__)
    return atomicOr(destination, value);
#elif defined(_MSC_VER)
    return uint32_t(_InterlockedOr(reinterpret_cast<volatile long *>(destination),
                                  static_cast<long>(value)));
#else
    return __atomic_fetch_or(destination, value, __ATOMIC_RELAXED);
#endif
}

YBI_CPU_GPU uint32_t AtomicAdd(uint32_t *destination, uint32_t value)
{
#if defined(__CUDA_ARCH__)
    return atomicAdd(destination, value);
#elif defined(_MSC_VER)
    return uint32_t(_InterlockedExchangeAdd(reinterpret_cast<volatile long *>(destination),
                                           static_cast<long>(value)));
#else
    return __atomic_fetch_add(destination, value, __ATOMIC_RELAXED);
#endif
}

#if defined(__CUDACC__)

YBI_GPU void AtomicMin(float *destination, float value)
{
    const int bits = __float_as_int(value);
    if (bits >= 0)
        atomicMin(reinterpret_cast<int *>(destination), bits);
    else
        atomicMax(reinterpret_cast<unsigned int *>(destination), __float_as_uint(value));
}

YBI_GPU void AtomicMax(float *destination, float value)
{
    const int bits = __float_as_int(value);
    if (bits >= 0)
        atomicMax(reinterpret_cast<int *>(destination), bits);
    else
        atomicMin(reinterpret_cast<unsigned int *>(destination), __float_as_uint(value));
}

#endif

} // namespace ybi
