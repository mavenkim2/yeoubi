#pragma once

#include "core/affine_transform.h"
#include "device/cuda_compat.h"
#include "device/device_texture.h"
#include "kernel/sampling/types.h"

namespace ybi
{

struct DeviceEnvmapLight
{
    AffineTransform worldFromLocal = AffineTransform::Identity();
    AffineTransform localFromWorld = AffineTransform::Identity();
    uint64_t textureObject = 0;
    TextureFormat textureFormat = TextureFormat::RGBA32_FLOAT;
    const AliasEntry *conditionalAlias = nullptr;
    const AliasEntry *marginalAlias = nullptr;
    uint32_t width = 0;
    uint32_t height = 0;
    uint32_t distributionWidth = 0;
    uint32_t distributionHeight = 0;
};

} // namespace ybi
