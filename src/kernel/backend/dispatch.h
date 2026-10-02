#pragma once

#include "device/cuda_compat.h"

#if defined(__CUDACC__) && defined(OPTIX_VERSION)
#include <optix_device.h>
#elif defined(__CUDACC__)
#error "DeviceDispatch requires OptiX; CUDA-only dispatch is not implemented."
#else
#error "DeviceDispatch requires OptiX device compilation; no supported device backend is active."
#endif

namespace ybi
{

template <typename Return, typename... Args>
YBI_GPU Return DeviceDispatch(uint32_t callableIndex, Args... args)
{
    return optixDirectCall<Return>(callableIndex, args...);
}

} // namespace ybi
