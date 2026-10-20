#pragma once

#include "core/affine_transform.h"
#include "device/cuda_compat.h"
#include "device/device_texture.h"

namespace ybi
{

struct DeviceEnvmapLight
{
    AffineTransform worldFromLocal = AffineTransform::Identity();
    AffineTransform localFromWorld = AffineTransform::Identity();
    uint64_t textureObject = 0;
    TextureFormat textureFormat = TextureFormat::RGBA32_FLOAT;
    const float *conditionalCDF = nullptr;
    const float *marginalCDF = nullptr;
    uint32_t width = 0;
    uint32_t height = 0;
};

} // namespace ybi
