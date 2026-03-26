#pragma once

#include "texture/virtual_texture/pixel_format.h"

#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

namespace ybi
{
namespace tilebin
{

struct TileFileHeader
{
    char magic[8];
    uint32_t version = 4;
    uint32_t channels = 4;
    uint32_t elementType = 1; // 1 = float32
    uint32_t udimCount = 0;
    uint64_t udimTableOffset = 0;
};

struct UdimEntry
{
    uint32_t udim = 0;
    uint32_t imageWidth = 0;
    uint32_t imageHeight = 0;
    uint32_t tileSize = 0;
    uint32_t mipCount = 0;
    uint32_t streamMipCount = 0;
    uint32_t tailMipCount = 0;
    uint32_t pixelFormat = static_cast<uint32_t>(texture::VirtualTexturePixelFormat::RGBA32_FLOAT);
    uint64_t mipRecordOffset = 0;
    uint32_t mipRecordCount = 0;
    uint32_t reserved1 = 0;
    uint64_t tileRecordOffset = 0;
    uint32_t tileRecordCount = 0;
    uint32_t reserved2 = 0;
    uint64_t payloadOffset = 0;
    uint64_t payloadBytes = 0;
};

struct MipRecord
{
    uint32_t mipLevel = 0;
    uint32_t width = 0;
    uint32_t height = 0;
    uint32_t tileCountX = 0;
    uint32_t tileCountY = 0;
    uint32_t firstTileRecord = 0;
    uint32_t tileRecordCount = 0;
    uint32_t isTail = 0;
    uint64_t byteOffset = 0;
    uint64_t byteSize = 0;
};

struct TileRecord
{
    uint32_t mipLevel = 0;
    uint32_t tileX = 0;
    uint32_t tileY = 0;
    uint32_t width = 0;
    uint32_t height = 0;
    uint64_t byteOffset = 0;
    uint64_t byteSize = 0;
};

struct UdimMipImage
{
    uint32_t level = 0;
    uint32_t width = 0;
    uint32_t height = 0;
    texture::VirtualTexturePixelFormat pixelFormat = texture::VirtualTexturePixelFormat::RGBA32_FLOAT;
    std::vector<float> rgba;
};

struct UdimImage
{
    uint32_t udim = 0;
    uint32_t width = 0;
    uint32_t height = 0;
    texture::VirtualTexturePixelFormat pixelFormat = texture::VirtualTexturePixelFormat::RGBA32_FLOAT;
    std::vector<UdimMipImage> mipLevels;
};

struct DiffStats
{
    float maxAbs = 0.0f;
    double meanAbs = 0.0;
    double rmse = 0.0;
    uint64_t mismatchCount = 0;
    size_t firstMismatch = 0;
};

bool WriteTileBinary(const std::filesystem::path &path,
                     int tileSize,
                     const std::vector<UdimImage> &images,
                     std::string *outError);

bool ReadTileBinary(const std::filesystem::path &path,
                    TileFileHeader &outHeader,
                    std::vector<UdimEntry> &outEntries,
                    std::vector<UdimImage> &outImages,
                    std::string *outError);

bool DiffImagesExact(const std::vector<float> &a,
                     const std::vector<float> &b,
                     float eps,
                     DiffStats *outStats);

} // namespace tilebin
} // namespace ybi
