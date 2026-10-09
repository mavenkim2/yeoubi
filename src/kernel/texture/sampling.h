#pragma once

#include "core/vec4f.h"
#include "device/device_texture.h"

#ifdef __CUDA_ARCH__
#include "kernel/backend/cuda/texture_sampling.h"
#else
#endif

namespace ybi
{

YBI_GPU float SampleDeviceTexture1D(uint64_t handle,
                                    uint32_t width,
                                    float u)
{
#ifdef __CUDA_ARCH__
    return SampleCUDATexture1D(handle, u);
#else
    const float *values = reinterpret_cast<const float *>(handle);
    const float position = Clamp(u * width - 0.5f, 0.f, float(width - 1));
    const uint32_t low = static_cast<uint32_t>(position);
    const uint32_t high = Min(low + 1, width - 1);
    return Lerp(values[low], values[high], position - low);
#endif
}

YBI_GPU Vec4f
SampleDeviceTexture(uint64_t handle, TextureFormat format, float u, float v, uint32_t layer)
{
#ifdef __CUDA_ARCH__
    return SampleCUDATexture(handle, format, u, v, layer);
#else
    return Vec4f(0.f);
#endif
}

YBI_GPU Vec4f
SampleDeviceTexture3D(uint64_t handle, TextureFormat format, float u, float v, float w)
{
#ifdef __CUDA_ARCH__
    return SampleCUDATexture3D(handle, format, u, v, w);
#else
    return Vec4f(0.f);
#endif
}

} // namespace ybi
