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
    uint32_t version = 1;
    uint32_t tileSize = 0;
    uint32_t imageWidth = 0;
    uint32_t imageHeight = 0;
    uint32_t channels = 4;
    uint32_t tileCountX = 0;
    uint32_t tileCountY = 0;
    uint32_t tileCount = 0;
    uint32_t elementType = 1; // 1 = float32
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

struct DiffStats
{
    float maxAbs = 0.0f;
    double meanAbs = 0.0;
    double rmse = 0.0;
    uint64_t mismatchCount = 0;
    size_t firstMismatch = 0;
};

bool WriteTileBinary(const std::filesystem::path &path,
                     int imageWidth,
                     int imageHeight,
                     int tileSize,
                     const std::vector<float> &image,
                     std::string *outError);

bool ReadTileBinary(const std::filesystem::path &path,
                    TileFileHeader &outHeader,
                    std::vector<TileRecord> &outRecords,
                    std::vector<float> &outImage,
                    std::string *outError);

bool DiffImagesExact(const std::vector<float> &a,
                     const std::vector<float> &b,
                     float eps,
                     DiffStats *outStats);

} // namespace tilebin
} // namespace ybi
