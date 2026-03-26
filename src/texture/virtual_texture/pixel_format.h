#pragma once

#include <cstddef>
#include <cstdint>

namespace ybi
{
namespace texture
{

enum class VirtualTexturePixelFormat : uint32_t
{
    RGBA32_FLOAT = 0,
    R16_FLOAT = 1,
    R32_FLOAT = 2,
    RG16_FLOAT = 3,
    RG32_FLOAT = 4,
    RGBA16_FLOAT = 5,
};

inline uint32_t VirtualTexturePixelFormatChannelCount(VirtualTexturePixelFormat format)
{
    switch (format)
    {
        case VirtualTexturePixelFormat::R16_FLOAT:
        case VirtualTexturePixelFormat::R32_FLOAT:
            return 1u;
        case VirtualTexturePixelFormat::RG16_FLOAT:
        case VirtualTexturePixelFormat::RG32_FLOAT:
            return 2u;
        case VirtualTexturePixelFormat::RGBA16_FLOAT:
        case VirtualTexturePixelFormat::RGBA32_FLOAT:
            return 4u;
        default:
            return 0u;
    }
}

inline uint32_t VirtualTexturePixelFormatBytesPerChannel(VirtualTexturePixelFormat format)
{
    switch (format)
    {
        case VirtualTexturePixelFormat::R16_FLOAT:
        case VirtualTexturePixelFormat::RG16_FLOAT:
        case VirtualTexturePixelFormat::RGBA16_FLOAT:
            return 2u;
        case VirtualTexturePixelFormat::R32_FLOAT:
        case VirtualTexturePixelFormat::RG32_FLOAT:
        case VirtualTexturePixelFormat::RGBA32_FLOAT:
            return 4u;
        default:
            return 0u;
    }
}

inline size_t VirtualTexturePixelFormatBytesPerPixel(VirtualTexturePixelFormat format)
{
    return static_cast<size_t>(VirtualTexturePixelFormatChannelCount(format)) *
           static_cast<size_t>(VirtualTexturePixelFormatBytesPerChannel(format));
}

inline bool IsValidVirtualTexturePixelFormat(const uint32_t value)
{
    switch (static_cast<VirtualTexturePixelFormat>(value))
    {
        case VirtualTexturePixelFormat::RGBA32_FLOAT:
        case VirtualTexturePixelFormat::R16_FLOAT:
        case VirtualTexturePixelFormat::R32_FLOAT:
        case VirtualTexturePixelFormat::RG16_FLOAT:
        case VirtualTexturePixelFormat::RG32_FLOAT:
        case VirtualTexturePixelFormat::RGBA16_FLOAT:
            return true;
        default:
            return false;
    }
}

} // namespace texture
} // namespace ybi
