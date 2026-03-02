#pragma once

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
    uint32_t version = 2;
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
    uint32_t tileCountX = 0;
    uint32_t tileCountY = 0;
    uint32_t tileCount = 0;
    uint64_t tileRecordOffset = 0;
    uint32_t tileRecordCount = 0;
    uint32_t reserved0 = 0;
    uint64_t payloadOffset = 0;
    uint64_t payloadBytes = 0;
};

struct TileRecord
{
    uint32_t tileX = 0;
    uint32_t tileY = 0;
    uint32_t width = 0;
    uint32_t height = 0;
    uint64_t byteOffset = 0;
    uint64_t byteSize = 0;
};

struct UdimImage
{
    uint32_t udim = 0;
    uint32_t width = 0;
    uint32_t height = 0;
    std::vector<float> rgba;
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
