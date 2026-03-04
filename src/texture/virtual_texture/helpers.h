#pragma once

#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace ybi
{
namespace texture
{

std::string ResolveVirtualTextureTileBinPath(const std::string &tilesDir,
                                             const std::string &texturePath);

std::unordered_map<unsigned int, std::string>
BuildVirtualTextureTileFileMap(const std::vector<std::pair<unsigned int, std::string>> &texturePaths,
                               const std::string &tilesDir);

} // namespace texture
} // namespace ybi
