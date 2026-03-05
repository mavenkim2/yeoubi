#pragma once

#include "tile_binary.h"

#include <string>
#include <vector>

namespace ybi
{
namespace tileprep
{

struct ExrMipChainResult
{
    bool hasStoredMipLevels = false;
    int tiled = 0;
    int tileLevelMode = -1;
    std::vector<tilebin::UdimMipImage> mipLevels;
};

bool LoadExrMipChain(const std::string &path,
                     bool flipVertical,
                     ExrMipChainResult *outResult,
                     std::string *outError);

} // namespace tileprep
} // namespace ybi
