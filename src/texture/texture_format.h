#pragma once

#include <cstddef>
#include <cstdint>

namespace ybi
{
enum class TextureFormat : uint32_t
{
    RGBA8_UNORM = 0,
    RGBA16_FLOAT = 1,
    R16_FLOAT = 2,
    R32_FLOAT = 3,
    RG16_FLOAT = 4,
    RG32_FLOAT = 5,
    RGBA32_FLOAT = 6,
};

inline uint32_t TextureFormatChannelCount(TextureFormat format)
{
    switch (format)
    {
        case TextureFormat::R16_FLOAT:
        case TextureFormat::R32_FLOAT:
            return 1u;
        case TextureFormat::RG16_FLOAT:
        case TextureFormat::RG32_FLOAT:
            return 2u;
        case TextureFormat::RGBA8_UNORM:
        case TextureFormat::RGBA16_FLOAT:
        case TextureFormat::RGBA32_FLOAT:
            return 4u;
        default:
            return 0u;
    }
}

inline uint32_t TextureFormatBytesPerChannel(TextureFormat format)
{
    switch (format)
    {
        case TextureFormat::RGBA8_UNORM:
            return 1u;
        case TextureFormat::R16_FLOAT:
        case TextureFormat::RG16_FLOAT:
        case TextureFormat::RGBA16_FLOAT:
            return 2u;
        case TextureFormat::R32_FLOAT:
        case TextureFormat::RG32_FLOAT:
        case TextureFormat::RGBA32_FLOAT:
            return 4u;
        default:
            return 0u;
    }
}

inline size_t TextureFormatPixelBytes(TextureFormat format)
{
    return static_cast<size_t>(TextureFormatChannelCount(format)) *
           static_cast<size_t>(TextureFormatBytesPerChannel(format));
}

inline bool IsValidTextureFormat(uint32_t value)
{
    switch (static_cast<TextureFormat>(value))
    {
        case TextureFormat::RGBA8_UNORM:
        case TextureFormat::RGBA16_FLOAT:
        case TextureFormat::R16_FLOAT:
        case TextureFormat::R32_FLOAT:
        case TextureFormat::RG16_FLOAT:
        case TextureFormat::RG32_FLOAT:
        case TextureFormat::RGBA32_FLOAT:
            return true;
        default:
            return false;
    }
}
} // namespace ybi
