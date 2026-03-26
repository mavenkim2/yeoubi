#pragma once

#include "tile_binary.h"

#include <string>
#include <vector>

namespace ybi
{
namespace tileprep
{

enum class ExrNumericType : uint8_t
{
    Unknown = 0,
    Float16 = 1,
    Float32 = 2,
};

struct ExrMipChainResult
{
    bool hasStoredMipLevels = false;
    int tiled = 0;
    int tileLevelMode = -1;
    uint32_t sourceChannelCount = 0u;
    ExrNumericType sourceNumericType = ExrNumericType::Unknown;
    std::vector<tilebin::UdimMipImage> mipLevels;
};

bool LoadExrMipChain(const std::string &path,
                     bool flipVertical,
                     ExrMipChainResult *outResult,
                     std::string *outError);

} // namespace tileprep
} // namespace ybi
