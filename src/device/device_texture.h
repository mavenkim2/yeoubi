#pragma once

#include "core/assert.h"
#include <cstdint>

namespace ybi
{

enum class TextureFormat
{
    RGBA8_UNORM,
    R16_FLOAT,
    RG16_FLOAT,
    RGBA16_FLOAT,
    R32_FLOAT,
    RG32_FLOAT,
    RGBA32_FLOAT,
};

enum class FilterType
{
    Point,
    Linear,
};

struct DeviceTextureInfo
{
    // TODO: 3D textures (not 2D texture arrays)
    TextureFormat format;
    FilterType filterType;
    uint32_t width = 0;
    uint32_t height = 0;
    uint32_t layers = 0;
    uint32_t mips = 0;
    // bool enableSurfaceWrites = true;
};

struct DeviceTexture
{
    uintptr_t array = 0;
    uint64_t textureObject = 0;
    DeviceTextureInfo info;
};

inline uint32_t GetTextureFormatSize(TextureFormat format)
{
    switch (format)
    {
        case TextureFormat::R16_FLOAT:
        {
            return 2;
        }
        case TextureFormat::RGBA8_UNORM:
        case TextureFormat::RG16_FLOAT:
        case TextureFormat::R32_FLOAT:
        {
            return 4;
        }
        case TextureFormat::RGBA16_FLOAT:
        case TextureFormat::RG32_FLOAT:
        {
            return 8;
        }
        case TextureFormat::RGBA32_FLOAT:
        {
            return 16;
        }
        default:
        {
            YBI_ASSERT(false);
            return 0;
        }
    }
}

} // namespace ybi
