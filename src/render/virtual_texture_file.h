#pragma once

#include "core/vec4f.h"
#include "device/device_texture.h"
#include <cstdint>

namespace ybi
{

#define YBVT_MAGIC 0x54564259u // YBVT

struct VirtualTextureFileHeader
{
    uint32_t magic;
    uint32_t version;

    uint32_t pageSize;
    uint32_t width;
    uint32_t height;
    uint32_t numMips;
    TextureFormat format;
    uint32_t compression;
    uint32_t numPages;
    Vec4f fallbackColor;

    uint64_t pageTableOffset;
    uint64_t pageDataOffset;
};

struct VirtualTexturePageEntry
{
    uint64_t offset;
    uint32_t packedSize;
    uint32_t unpackedSize;
};

} // namespace ybi
