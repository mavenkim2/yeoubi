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

void ExtractTileFloatSamples(const std::vector<float> &image,
                             int imageWidth,
                             int imageHeight,
                             uint32_t channelCount,
                             int tileX,
                             int tileY,
                             int tileSize,
                             std::vector<float> &outTile,
                             int &outWidth,
                             int &outHeight);

bool ValidateImageBasics(const std::vector<UdimImage> &images, int tileSize, std::string *outError);

bool ReadBytes(std::ifstream *stream, void *dst, size_t bytes);

} // namespace detail
} // namespace tilebin
} // namespace ybi
