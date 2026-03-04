#include "texture/virtual_texture/helpers.h"

#include "texture/udim_utils.h"

#include <cctype>
#include <filesystem>

namespace ybi
{
namespace texture
{
namespace
{

std::string SanitizeTileStem(const std::string &path)
{
    std::string out = path;
    for (char &c : out)
    {
        if (!(std::isalnum(static_cast<unsigned char>(c)) || c == '_' || c == '-' || c == '.'))
        {
            c = '_';
        }
    }
    if (out.empty())
    {
        out = "texture";
    }
    return out;
}

} // namespace

std::string ResolveVirtualTextureTileBinPath(const std::string &tilesDir,
                                             const std::string &texturePath)
{
    const std::string baseNoUdim = ybi::usd_ntc::StripUdimFromPath(texturePath);
    const std::string stem = SanitizeTileStem(baseNoUdim);
    return (std::filesystem::path(tilesDir) / (stem + ".tiles.bin")).string();
}

std::unordered_map<unsigned int, std::string>
BuildVirtualTextureTileFileMap(const std::vector<std::pair<unsigned int, std::string>> &texturePaths,
                               const std::string &tilesDir)
{
    std::unordered_map<unsigned int, std::string> out;
    if (tilesDir.empty())
    {
        return out;
    }

    out.reserve(texturePaths.size());
    for (const auto &entry : texturePaths)
    {
        const unsigned int textureId = entry.first;
        const std::string &texturePath = entry.second;
        if (texturePath.empty())
        {
            continue;
        }

        out.emplace(textureId, ResolveVirtualTextureTileBinPath(tilesDir, texturePath));
    }

    return out;
}

} // namespace texture
} // namespace ybi
