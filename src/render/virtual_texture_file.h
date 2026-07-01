#pragma once

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
    uint32_t numChannels;
    uint32_t numBytesPerChannel;
    uint32_t compression;
    uint32_t numPages;

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
