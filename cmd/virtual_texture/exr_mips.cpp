#include "exr_mips.h"
#include "texture/exr_io.h"
#include "third_party/tinyexr/tinyexr.h"
#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <map>
#include <string>
#include <vector>
namespace ybi
{
namespace tileprep
{
namespace
{
struct LayerChannels
{
    int first = -1;
    int r = -1;
    int g = -1;
    int b = -1;
    int a = -1;
};
std::pair<std::string, std::string> SplitChannelName(const char *name)
{
    const std::string full = name ? std::string(name) : std::string();
    const size_t dot = full.rfind('.');
    if (dot == std::string::npos)
    {
        return {"", full};
    }
    return {full.substr(0, dot), full.substr(dot + 1)};
}
void FlipVerticalRgba(std::vector<float> *rgba, uint32_t width, uint32_t height)
{
    if (!rgba || width == 0u || height <= 1u)
    {
        return;
    }
    const size_t rowFloats = static_cast<size_t>(width) * 4u;
    std::vector<float> tmp(rowFloats);
    for (uint32_t y = 0u; y < height / 2u; ++y)
    {
        float *rowA = rgba->data() + static_cast<size_t>(y) * rowFloats;
        float *rowB = rgba->data() + static_cast<size_t>(height - 1u - y) * rowFloats;
        std::memcpy(tmp.data(), rowA, rowFloats * sizeof(float));
        std::memcpy(rowA, rowB, rowFloats * sizeof(float));
        std::memcpy(rowB, tmp.data(), rowFloats * sizeof(float));
    }
}
void Downsample2x2(const std::vector<float> &src,
                   uint32_t srcW,
                   uint32_t srcH,
                   std::vector<float> *dst,
                   uint32_t *dstW,
                   uint32_t *dstH)
{
    *dstW = std::max(1u, srcW >> 1u);
    *dstH = std::max(1u, srcH >> 1u);
    dst->assign(static_cast<size_t>(*dstW) * static_cast<size_t>(*dstH) * 4u, 0.0f);
    for (uint32_t y = 0; y < *dstH; ++y)
    {
        for (uint32_t x = 0; x < *dstW; ++x)
        {
            const uint32_t sx0 = std::min(srcW - 1u, x * 2u + 0u);
            const uint32_t sx1 = std::min(srcW - 1u, x * 2u + 1u);
            const uint32_t sy0 = std::min(srcH - 1u, y * 2u + 0u);
            const uint32_t sy1 = std::min(srcH - 1u, y * 2u + 1u);
            const size_t outBase =
                (static_cast<size_t>(y) * static_cast<size_t>(*dstW) + static_cast<size_t>(x)) *
                4u;
            for (uint32_t c = 0; c < 4u; ++c)
            {
                const size_t i00 = (static_cast<size_t>(sy0) * static_cast<size_t>(srcW) +
                                    static_cast<size_t>(sx0)) *
                                       4u +
                                   c;
                const size_t i10 = (static_cast<size_t>(sy0) * static_cast<size_t>(srcW) +
                                    static_cast<size_t>(sx1)) *
                                       4u +
                                   c;
                const size_t i01 = (static_cast<size_t>(sy1) * static_cast<size_t>(srcW) +
                                    static_cast<size_t>(sx0)) *
                                       4u +
                                   c;
                const size_t i11 = (static_cast<size_t>(sy1) * static_cast<size_t>(srcW) +
                                    static_cast<size_t>(sx1)) *
                                       4u +
                                   c;
                (*dst)[outBase + c] = 0.25f * (src[i00] + src[i10] + src[i01] + src[i11]);
            }
        }
    }
}
void BuildManualMipChain(const std::vector<float> &base,
                         uint32_t width,
                         uint32_t height,
                         std::vector<tilebin::UdimMipImage> *outMips)
{
    outMips->clear();
    tilebin::UdimMipImage mip0 = {};
    mip0.level = 0u;
    mip0.width = width;
    mip0.height = height;
    mip0.rgba = base;
    outMips->push_back(std::move(mip0));
    while (outMips->back().width > 1u || outMips->back().height > 1u)
    {
        tilebin::UdimMipImage next = {};
        next.level = static_cast<uint32_t>(outMips->size());
        Downsample2x2(outMips->back().rgba,
                      outMips->back().width,
                      outMips->back().height,
                      &next.rgba,
                      &next.width,
                      &next.height);
        outMips->push_back(std::move(next));
    }
}
bool FindRgbaChannels(
    const EXRHeader &header, int *outR, int *outG, int *outB, int *outA, std::string *outError)
{
    std::map<std::string, LayerChannels> layers;
    for (int i = 0; i < header.num_channels; ++i)
    {
        const auto split = SplitChannelName(header.channels[i].name);
        printf("%s name: %s index: %i\n", header.name, header.channels[i].name, i);
        LayerChannels &layer = layers[split.first];
        if (layer.first < 0)
        {
            layer.first = i;
        }
        if (split.second == "R" && layer.r < 0)
        {
            layer.r = i;
        }
        else if (split.second == "G" && layer.g < 0)
        {
            layer.g = i;
        }
        else if (split.second == "B" && layer.b < 0)
        {
            layer.b = i;
        }
        else if (split.second == "A" && layer.a < 0)
        {
            layer.a = i;
        }
    }
    auto choose = layers.find("");
    if (choose == layers.end() || choose->second.r < 0 || choose->second.g < 0 ||
        choose->second.b < 0)
    {
        choose = layers.end();
        for (auto it = layers.begin(); it != layers.end(); ++it)
        {
            if (it->second.r >= 0 && it->second.g >= 0 && it->second.b >= 0)
            {
                choose = it;
                break;
            }
        }
    }
    if (choose == layers.end())
    {
        choose = layers.begin();
        for (auto it = layers.begin(); it != layers.end(); ++it)
        {
            if (it->second.first >= 0)
            {
                choose = it;
                break;
            }
        }
    }
    if (choose == layers.end() || choose->second.first < 0)
    {
        if (outError)
        {
            *outError = "could not find image channels";
        }
        return false;
    }
    *outR = (choose->second.r >= 0) ? choose->second.r : choose->second.first;
    *outG = (choose->second.g >= 0) ? choose->second.g : choose->second.first;
    *outB = (choose->second.b >= 0) ? choose->second.b : choose->second.first;
    *outA = choose->second.a;
    return true;
}
bool ReadLevelRgba(const EXRHeader &header,
                   const EXRImage &level,
                   int idxR,
                   int idxG,
                   int idxB,
                   int idxA,
                   std::vector<float> *outRgba,
                   std::string *outError)
{
    if (level.width <= 0 || level.height <= 0)
    {
        if (outError)
        {
            *outError = "invalid EXR level dimensions";
        }
        return false;
    }
    const uint32_t width = static_cast<uint32_t>(level.width);
    const uint32_t height = static_cast<uint32_t>(level.height);
    outRgba->assign(static_cast<size_t>(width) * static_cast<size_t>(height) * 4u, 1.0f);
    if (level.tiles && level.num_tiles > 0)
    {
        const int tileSizeX = std::max(1, header.tile_size_x);
        const int tileSizeY = std::max(1, header.tile_size_y);
        for (int t = 0; t < level.num_tiles; ++t)
        {
            const EXRTile &tile = level.tiles[t];
            if (!tile.images)
            {
                continue;
            }
            float **channels = reinterpret_cast<float **>(tile.images);
            if (!channels[idxR] || !channels[idxG] || !channels[idxB])
            {
                continue;
            }
            const int tw = std::max(1, tile.width > 0 ? tile.width : tileSizeX);
            const int th = std::max(1, tile.height > 0 ? tile.height : tileSizeY);
            const int x0 = tile.offset_x * tileSizeX;
            const int y0 = tile.offset_y * tileSizeY;
            for (int y = 0; y < th; ++y)
            {
                const int dstY = y0 + y;
                if (dstY < 0 || dstY >= level.height)
                {
                    continue;
                }
                for (int x = 0; x < tw; ++x)
                {
                    const int dstX = x0 + x;
                    if (dstX < 0 || dstX >= level.width)
                    {
                        continue;
                    }
                    // TinyEXR stores each tile row at the header tile stride, even
                    // when edge tiles are clipped smaller at this mip level.
                    const size_t srcIndex =
                        static_cast<size_t>(y) * static_cast<size_t>(tileSizeX) +
                        static_cast<size_t>(x);
                    const size_t dstBase =
                        (static_cast<size_t>(dstY) * static_cast<size_t>(width) +
                         static_cast<size_t>(dstX)) *
                        4u;
                    (*outRgba)[dstBase + 0u] = channels[idxR][srcIndex];
                    (*outRgba)[dstBase + 1u] = channels[idxG][srcIndex];
                    (*outRgba)[dstBase + 2u] = channels[idxB][srcIndex];
                    (*outRgba)[dstBase + 3u] =
                        (idxA >= 0 && channels[idxA]) ? channels[idxA][srcIndex] : 1.0f;
                }
            }
        }
        return true;
    }
    if (!level.images)
    {
        if (outError)
        {
            *outError = "EXR level has neither tiles nor scanline image data";
        }
        return false;
    }
    float **channels = reinterpret_cast<float **>(level.images);
    if (!channels[idxR] || !channels[idxG] || !channels[idxB])
    {
        if (outError)
        {
            *outError = "EXR RGB channels are null";
        }
        return false;
    }
    const size_t pixelCount = static_cast<size_t>(width) * static_cast<size_t>(height);
    for (size_t i = 0; i < pixelCount; ++i)
    {
        const size_t dstBase = i * 4u;
        (*outRgba)[dstBase + 0u] = channels[idxR][i];
        (*outRgba)[dstBase + 1u] = channels[idxG][i];
        (*outRgba)[dstBase + 2u] = channels[idxB][i];
        (*outRgba)[dstBase + 3u] = (idxA >= 0 && channels[idxA]) ? channels[idxA][i] : 1.0f;
    }
    return true;
}
bool LoadStoredMipChain(const std::string &path,
                        bool flipVertical,
                        EXRHeader *header,
                        std::vector<tilebin::UdimMipImage> *outMips,
                        std::string *outError)
{
    for (int i = 0; i < header->num_channels; ++i)
    {
        if (header->pixel_types[i] == TINYEXR_PIXELTYPE_HALF ||
            header->pixel_types[i] == TINYEXR_PIXELTYPE_FLOAT)
        {
            header->requested_pixel_types[i] = TINYEXR_PIXELTYPE_FLOAT;
        }
    }
    EXRImage image = {};
    InitEXRImage(&image);
    const char *err = nullptr;
    const int loadRc = LoadEXRImageFromFile(&image, header, path.c_str(), &err);
    if (loadRc != TINYEXR_SUCCESS)
    {
        if (outError)
        {
            *outError = "LoadEXRImageFromFile failed";
            if (err)
            {
                *outError += " (";
                *outError += err;
                *outError += ")";
            }
        }
        if (err)
        {
            FreeEXRErrorMessage(err);
        }
        FreeEXRImage(&image);
        return false;
    }
    int idxR = -1;
    int idxG = -1;
    int idxB = -1;
    int idxA = -1;

    if (!FindRgbaChannels(*header, &idxR, &idxG, &idxB, &idxA, outError))
    {
        FreeEXRImage(&image);
        return false;
    }
    if (header->requested_pixel_types[idxR] != TINYEXR_PIXELTYPE_FLOAT ||
        header->requested_pixel_types[idxG] != TINYEXR_PIXELTYPE_FLOAT ||
        header->requested_pixel_types[idxB] != TINYEXR_PIXELTYPE_FLOAT ||
        (idxA >= 0 && header->requested_pixel_types[idxA] != TINYEXR_PIXELTYPE_FLOAT))
    {
        if (outError)
        {
            *outError = "EXR RGB(A) channels are not float-convertible";
        }
        FreeEXRImage(&image);
        return false;
    }
    std::vector<tilebin::UdimMipImage> loaded;
    const EXRImage *level = &image;
    while (level)
    {
        tilebin::UdimMipImage mip = {};
        mip.level = static_cast<uint32_t>(std::max(level->level_x, 0));
        mip.width = static_cast<uint32_t>(std::max(level->width, 1));
        mip.height = static_cast<uint32_t>(std::max(level->height, 1));
        if (!ReadLevelRgba(*header, *level, idxR, idxG, idxB, idxA, &mip.rgba, outError))
        {
            FreeEXRImage(&image);
            return false;
        }
        if (flipVertical)
        {
            FlipVerticalRgba(&mip.rgba, mip.width, mip.height);
        }
        loaded.push_back(std::move(mip));
        level = level->next_level;
    }
    FreeEXRImage(&image);
    if (loaded.empty())
    {
        if (outError)
        {
            *outError = "EXR has no mip levels";
        }
        return false;
    }
    if (loaded.front().level != 0u)
    {
        if (outError)
        {
            *outError = "EXR mip chain missing level 0";
        }
        return false;
    }
    outMips->clear();
    outMips->reserve(loaded.size());
    outMips->push_back(std::move(loaded[0]));
    size_t nextLoadedIndex = 1u;
    while (outMips->back().width > 1u || outMips->back().height > 1u)
    {
        const uint32_t nextLevel = static_cast<uint32_t>(outMips->size());
        const uint32_t expectedW = std::max(1u, outMips->back().width >> 1u);
        const uint32_t expectedH = std::max(1u, outMips->back().height >> 1u);

        // TinyEXR's mipmap loader appends next_level in ascending level order.
        while (nextLoadedIndex < loaded.size() && loaded[nextLoadedIndex].level < nextLevel)
        {
            ++nextLoadedIndex;
        }
        if (nextLoadedIndex < loaded.size() && loaded[nextLoadedIndex].level == nextLevel &&
            loaded[nextLoadedIndex].width == expectedW &&
            loaded[nextLoadedIndex].height == expectedH)
        {
            outMips->push_back(std::move(loaded[nextLoadedIndex]));
            ++nextLoadedIndex;
            continue;
        }
        tilebin::UdimMipImage next = {};
        next.level = nextLevel;
        Downsample2x2(outMips->back().rgba,
                      outMips->back().width,
                      outMips->back().height,
                      &next.rgba,
                      &next.width,
                      &next.height);
        outMips->push_back(std::move(next));
    }
    return true;
}
} // namespace
bool LoadExrMipChain(const std::string &path,
                     bool flipVertical,
                     ExrMipChainResult *outResult,
                     std::string *outError)
{
    if (!outResult)
    {
        if (outError)
        {
            *outError = "invalid output pointer";
        }
        return false;
    }
    *outResult = {};
    EXRVersion version = {};
    const int versionRc = ParseEXRVersionFromFile(&version, path.c_str());
    if (versionRc != TINYEXR_SUCCESS)
    {
        if (outError)
        {
            *outError = "ParseEXRVersionFromFile failed";
        }
        return false;
    }
    EXRHeader header = {};
    InitEXRHeader(&header);
    const char *err = nullptr;
    const int parseRc = ParseEXRHeaderFromFile(&header, &version, path.c_str(), &err);
    if (parseRc != TINYEXR_SUCCESS)
    {
        if (outError)
        {
            *outError = "ParseEXRHeaderFromFile failed";
            if (err)
            {
                *outError += " (";
                *outError += err;
                *outError += ")";
            }
        }
        if (err)
        {
            FreeEXRErrorMessage(err);
        }
        FreeEXRHeader(&header);
        return false;
    }
    outResult->tiled = header.tiled;
    outResult->tileLevelMode = header.tile_level_mode;
    outResult->hasStoredMipLevels =
        (header.tiled != 0 && header.tile_level_mode == TINYEXR_TILE_MIPMAP_LEVELS);
    bool ok = true;
    if (outResult->hasStoredMipLevels)
    {
        ok = LoadStoredMipChain(path, flipVertical, &header, &outResult->mipLevels, outError);
    }
    else
    {
        int width = 0;
        int height = 0;
        std::vector<float> base;
        std::string reason;
        if (!ybi::texture::LoadExrRgba(path, &width, &height, &base, &reason, flipVertical))
        {
            if (outError)
            {
                *outError = reason;
            }
            ok = false;
        }
        else
        {
            BuildManualMipChain(base,
                                static_cast<uint32_t>(width),
                                static_cast<uint32_t>(height),
                                &outResult->mipLevels);
        }
    }
    FreeEXRHeader(&header);
    return ok;
}
} // namespace tileprep
} // namespace ybi
