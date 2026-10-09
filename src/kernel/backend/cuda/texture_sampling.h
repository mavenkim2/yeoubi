#pragma once

#include "core/vec4f.h"
#include "device/device_texture.h"
#include <cuda_runtime.h>

namespace ybi
{

YBI_GPU float SampleCUDATexture1D(uint64_t handle, float u)
{
    const cudaTextureObject_t texture = static_cast<cudaTextureObject_t>(handle);
    return tex1D<float>(texture, u);
}

YBI_GPU Vec4f SampleCUDATexture(
    uint64_t handle,
    TextureFormat format,
    float u,
    float v,
    uint32_t layer)
{
    const cudaTextureObject_t texture =
        reinterpret_cast<cudaTextureObject_t>(handle);

    switch (GetTextureFormatNumChannels(format))
    {
        case 1:
        {
            const float x = tex2DLayered<float>(texture, u, v, layer);

            return Vec4f(x, x, x, 1.0f);
        }

        case 2:
        {
            const float2 xy =
                tex2DLayered<float2>(texture, u, v, layer);

            return Vec4f(xy.x, xy.y, 0.0f, 1.0f);
        }

        case 4:
        {
            const float4 xyzw =
                tex2DLayered<float4>(texture, u, v, layer);

            return Vec4f(xyzw.x, xyzw.y, xyzw.z, xyzw.w);
        }
    }

    YBI_ASSERT(false);
    return {};
}

YBI_GPU Vec4f SampleCUDATexture3D(
    uint64_t handle, TextureFormat format, float u, float v, float w)
{
    const cudaTextureObject_t texture = reinterpret_cast<cudaTextureObject_t>(handle);

    switch (GetTextureFormatNumChannels(format))
    {
        case 1:
        {
            const float x = tex3D<float>(texture, u, v, w);
            return Vec4f(x, x, x, 1.f);
        }
        case 2:
        {
            const float2 xy = tex3D<float2>(texture, u, v, w);
            return Vec4f(xy.x, xy.y, 0.f, 1.f);
        }
        case 4:
        {
            const float4 xyzw = tex3D<float4>(texture, u, v, w);
            return Vec4f(xyzw.x, xyzw.y, xyzw.z, xyzw.w);
        }
    }

    YBI_ASSERT(false);
    return {};
}

}
