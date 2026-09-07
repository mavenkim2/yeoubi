#pragma once

#include "core/vec4f.h"
#include "device/device_texture.h"

#ifdef __CUDA_ARCH__
#include "render/cuda/texture_sampling.h"
#else
#endif

namespace ybi
{

YBI_GPU Vec4f
SampleDeviceTexture(uint64_t handle, TextureFormat format, float u, float v, uint32_t layer)
{
#ifdef __CUDA_ARCH__
    SampleCUDATexture(handle, format, u, v, layer);
#else
#error Device not supported yet
#endif
}

} // namespace ybi
