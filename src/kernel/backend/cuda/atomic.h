#pragma once

#include "device/cuda_compat.h"

#if defined(__CUDACC__)

namespace ybi
{

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

} // namespace ybi

#endif
