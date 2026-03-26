#include "exr_mips.h"
#include "texture/exr_io.h"
#include "third_party/tinyexr/tinyexr.h"
#include "util/assert.h"
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
    std::vector<int> ordered;
    int r = -1;
    int g = -1;
    int b = -1;
    int a = -1;
};

struct SelectedChannels
{
    std::array<int, 4> rgba = {-1, -1, -1, -1};
    std::array<int, 4> source = {-1, -1, -1, -1};
    uint32_t sourceChannelCount = 0u;
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
void FlipImageVertical(std::vector<float> *values,
                       uint32_t width,
                       uint32_t height,
                       int numChannels)
{
    if (!values || width == 0u || height <= 1u)
    {
        return;
    }
    const size_t rowFloats = static_cast<size_t>(width) * numChannels;
    std::vector<float> tmp(rowFloats);
    for (uint32_t y = 0u; y < height / 2u; ++y)
    {
        float *rowA = values->data() + static_cast<size_t>(y) * rowFloats;
        float *rowB = values->data() + static_cast<size_t>(height - 1u - y) * rowFloats;
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
bool FindSelectedChannels(const EXRHeader &header,
                          SelectedChannels *outChannels,
                          std::string *outError)
{
    if (!outChannels)
    {
        return false;
    }
    *outChannels = {};
    std::map<std::string, LayerChannels> layers;
    for (int i = 0; i < header.num_channels; ++i)
    {
        const auto split = SplitChannelName(header.channels[i].name);
        LayerChannels &layer = layers[split.first];
        layer.ordered.push_back(i);
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
            if (!it->second.ordered.empty())
            {
                choose = it;
                break;
            }
        }
    }
    if (choose == layers.end() || choose->second.ordered.empty())
    {
        if (outError)
        {
            *outError = "could not find image channels";
        }
        return false;
    }

    if (choose->second.ordered.size() > outChannels->source.size())
    {
        if (outError)
        {
            *outError = "selected EXR layer has more than 4 channels";
        }
        return false;
    }

    outChannels->sourceChannelCount = static_cast<uint32_t>(choose->second.ordered.size());
    for (size_t i = 0; i < choose->second.ordered.size(); ++i)
    {
        outChannels->source[i] = choose->second.ordered[i];
    }

    const int first = choose->second.ordered[0];
    outChannels->rgba[0] = (choose->second.r >= 0) ? choose->second.r : first;
    outChannels->rgba[1] = (choose->second.g >= 0) ? choose->second.g : first;
    outChannels->rgba[2] = (choose->second.b >= 0) ? choose->second.b : first;
    outChannels->rgba[3] = choose->second.a;
    return true;
}

bool DetermineSourceNumericType(const EXRHeader &header,
                                const SelectedChannels &channels,
                                ExrNumericType *outType,
                                std::string *outError)
{
    if (!outType)
    {
        return false;
    }
    *outType = ExrNumericType::Unknown;
    if (channels.sourceChannelCount == 0u)
    {
        if (outError)
        {
            *outError = "selected EXR layer has no channels";
        }
        return false;
    }

    int pixelType = -1;
    for (uint32_t i = 0u; i < channels.sourceChannelCount; ++i)
    {
        const int channelIndex = channels.source[i];
        if (channelIndex < 0 || channelIndex >= header.num_channels)
        {
            if (outError)
            {
                *outError = "selected EXR channel index out of range";
            }
            return false;
        }

        const int currentType = header.pixel_types[channelIndex];
        if (currentType != TINYEXR_PIXELTYPE_HALF && currentType != TINYEXR_PIXELTYPE_FLOAT)
        {
            if (outError)
            {
                *outError = "selected EXR channels are not half/float";
            }
            return false;
        }
        if (pixelType < 0)
        {
            pixelType = currentType;
            continue;
        }
        if (pixelType != currentType)
        {
            if (outError)
            {
                *outError = "selected EXR channels use mixed numeric types";
            }
            return false;
        }
    }

    *outType = (pixelType == TINYEXR_PIXELTYPE_HALF) ? ExrNumericType::Float16
                                                     : ExrNumericType::Float32;
    return true;
}

bool ReadLevelRgba(const EXRHeader &header,
                   const EXRImage &level,
                   const std::array<int, 4> &channelIndices,
                   std::vector<float> *outValues,
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
    outValues->assign(static_cast<size_t>(width) * static_cast<size_t>(height) * 4u, 1.0f);
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
            for (int channel = 0; channel < 4; ++channel)
            {
                const int channelIndex = channelIndices[channel];
                if (channelIndex < 0)
                {
                    continue;
                }
                if (channelIndex >= header.num_channels || !channels[channelIndex])
                {
                    if (outError)
                    {
                        *outError = "EXR RGB channels are null";
                    }
                    return false;
                }
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
                    for (int channel = 0; channel < 4; ++channel)
                    {
                        const int channelIndex = channelIndices[channel];
                        if (channelIndex >= 0)
                        {
                            (*outValues)[dstBase + channel] = channels[channelIndex][srcIndex];
                        }
                        else if (channel == 3)
                        {
                            (*outValues)[dstBase + channel] = 1.0f;
                        }
                    }
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
    for (int channel = 0; channel < 4; ++channel)
    {
        const int channelIndex = channelIndices[channel];
        if (channelIndex < 0)
        {
            continue;
        }
        if (channelIndex >= header.num_channels || !channels[channelIndex])
        {
            if (outError)
            {
                *outError = "EXR RGB channels are null";
            }
            return false;
        }
    }
    const size_t pixelCount = static_cast<size_t>(width) * static_cast<size_t>(height);
    for (size_t i = 0; i < pixelCount; ++i)
    {
        const size_t dstBase = i * 4u;
        for (int channel = 0; channel < 4; ++channel)
        {
            const int channelIndex = channelIndices[channel];
            if (channelIndex >= 0)
            {
                (*outValues)[dstBase + channel] = channels[channelIndex][i];
            }
            else if (channel == 3)
            {
                (*outValues)[dstBase + channel] = 1.0f;
            }
        }
    }
    return true;
}
bool LoadStoredMipChain(const std::string &path,
                        bool flipVertical,
                        EXRHeader *header,
                        std::vector<tilebin::UdimMipImage> *outMips,
                        std::string *outError)
{
    YBI_ERROR(header->num_channels <= 4, "multi layer exrs not supported yet");
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
    SelectedChannels channels;

    if (!FindSelectedChannels(*header, &channels, outError))
    {
        FreeEXRImage(&image);
        return false;
    }
    for (int channel = 0; channel < header->num_channels; channel++)
    {
        if (header->requested_pixel_types[channel] != TINYEXR_PIXELTYPE_FLOAT)
        {
            if (outError)
            {
                *outError = "EXR RGB(A) channels are not float-convertible";
            }
            FreeEXRImage(&image);
            return false;
        }
    }
    std::vector<tilebin::UdimMipImage> loaded;
    const EXRImage *level = &image;
    while (level)
    {
        tilebin::UdimMipImage mip = {};
        mip.level = static_cast<uint32_t>(std::max(level->level_x, 0));
        mip.width = static_cast<uint32_t>(std::max(level->width, 1));
        mip.height = static_cast<uint32_t>(std::max(level->height, 1));
        if (!ReadLevelRgba(*header, *level, channels.rgba, &mip.rgba, outError))
        {
            FreeEXRImage(&image);
            return false;
        }
        if (flipVertical)
        {
            FlipImageVertical(&mip.rgba, mip.width, mip.height, 4);
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
    SelectedChannels sourceChannels = {};
    if (!FindSelectedChannels(header, &sourceChannels, outError))
    {
        FreeEXRHeader(&header);
        return false;
    }
    outResult->sourceChannelCount = sourceChannels.sourceChannelCount;
    if (!DetermineSourceNumericType(header, sourceChannels, &outResult->sourceNumericType, outError))
    {
        FreeEXRHeader(&header);
        return false;
    }
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
