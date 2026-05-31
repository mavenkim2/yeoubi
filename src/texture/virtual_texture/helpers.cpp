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
    const std::string baseNoUdim = ybi::texture::StripUdimFromPath(texturePath);
    const std::string suffix =
        SanitizeTileStem(std::filesystem::path(baseNoUdim).filename().string()) + ".tiles.bin";
    std::filesystem::path bestMatch = {};
    std::error_code ec;
    for (const std::filesystem::directory_entry &entry :
         std::filesystem::directory_iterator(std::filesystem::path(tilesDir), ec))
    {
        if (ec || !entry.is_regular_file())
        {
            continue;
        }
        const std::string filename = entry.path().filename().string();
        if (filename == suffix ||
            (filename.size() > suffix.size() &&
             filename.compare(filename.size() - suffix.size(), suffix.size(), suffix) == 0))
        {
            if (bestMatch.empty())
            {
                bestMatch = entry.path();
                continue;
            }
            const std::string bestName = bestMatch.filename().string();
            if (filename.size() < bestName.size() || (filename.size() == bestName.size() && filename < bestName))
            {
                bestMatch = entry.path();
            }
        }
    }
    if (!bestMatch.empty())
    {
        return bestMatch.string();
    }
    return {};
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
