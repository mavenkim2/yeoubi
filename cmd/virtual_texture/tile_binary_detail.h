#pragma once

#include "tile_binary.h"

#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

namespace ybi
{
namespace tilebin
{
namespace detail
{

struct V2UdimEntry
{
    uint32_t udim = 0;
    uint32_t imageWidth = 0;
    uint32_t imageHeight = 0;
    uint32_t tileSize = 0;
    uint32_t tileCountX = 0;
    uint32_t tileCountY = 0;
    uint32_t tileCount = 0;
    uint64_t tileRecordOffset = 0;
    uint32_t tileRecordCount = 0;
    uint32_t reserved0 = 0;
    uint64_t payloadOffset = 0;
    uint64_t payloadBytes = 0;
};

struct V2TileRecord
{
    uint32_t tileX = 0;
    uint32_t tileY = 0;
    uint32_t width = 0;
    uint32_t height = 0;
    uint64_t byteOffset = 0;
    uint64_t byteSize = 0;
};

void ExtractTileRgbaF32(const std::vector<float> &image,
                        int imageWidth,
                        int imageHeight,
                        int tileX,
                        int tileY,
                        int tileSize,
                        std::vector<float> &outTile,
                        int &outWidth,
                        int &outHeight);

bool ValidateImageBasics(const std::vector<UdimImage> &images, int tileSize, std::string *outError);

bool ReadBytes(std::ifstream *stream, void *dst, size_t bytes);

bool ReadTileBinaryV2(std::ifstream *in,
                      const std::filesystem::path &path,
                      TileFileHeader &outHeader,
                      std::vector<UdimEntry> &outEntries,
                      std::vector<UdimImage> &outImages,
                      std::string *outError);

} // namespace detail
} // namespace tilebin
} // namespace ybi
