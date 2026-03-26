#include "tile_binary_detail.h"

#include <algorithm>
#include <cstring>

namespace ybi
{
namespace tilebin
{
namespace detail
{

namespace
{

static constexpr uint32_t kUdimMin = 1001;
static constexpr uint32_t kUdimMax = 1100;

} // namespace

void ExtractTileFloatSamples(const std::vector<float> &image,
                             int imageWidth,
                             int imageHeight,
                             uint32_t channelCount,
                             int tileX,
                             int tileY,
                             int tileSize,
                             std::vector<float> &outTile,
                             int &outWidth,
                             int &outHeight)
{
    const int x0 = tileX * tileSize;
    const int y0 = tileY * tileSize;
    outWidth = std::max(0, std::min(tileSize, imageWidth - x0));
    outHeight = std::max(0, std::min(tileSize, imageHeight - y0));
    outTile.assign(static_cast<size_t>(outWidth) * static_cast<size_t>(outHeight) *
                       static_cast<size_t>(channelCount),
                   0.0f);
    for (int y = 0; y < outHeight; ++y)
    {
        const float *src = image.data() + (static_cast<size_t>(y0 + y) * static_cast<size_t>(imageWidth) +
                                           static_cast<size_t>(x0)) *
                                              static_cast<size_t>(channelCount);
        float *dst = outTile.data() + static_cast<size_t>(y) * static_cast<size_t>(outWidth) *
                                          static_cast<size_t>(channelCount);
        std::memcpy(dst,
                    src,
                    static_cast<size_t>(outWidth) * static_cast<size_t>(channelCount) *
                        sizeof(float));
    }
}

bool ValidateImageBasics(const std::vector<UdimImage> &images, int tileSize, std::string *outError)
{
    if (tileSize <= 0)
    {
        if (outError)
        {
            *outError = "tileSize must be > 0";
        }
        return false;
    }
    if (images.empty())
    {
        if (outError)
        {
            *outError = "no UDIM images to write";
        }
        return false;
    }

    uint32_t prevUdim = 0u;
    bool havePrevUdim = false;
    for (const UdimImage &img : images)
    {
        if (img.udim < kUdimMin || img.udim > kUdimMax)
        {
            if (outError)
            {
                *outError = "invalid UDIM id";
            }
            return false;
        }
        if (havePrevUdim && img.udim <= prevUdim)
        {
            if (outError)
            {
                *outError = "UDIM ids must be strictly ascending";
            }
            return false;
        }
        prevUdim = img.udim;
        havePrevUdim = true;
        if (img.width == 0 || img.height == 0)
        {
            if (outError)
            {
                *outError = "invalid image dimensions";
            }
            return false;
        }
        const size_t bytesPerPixel = TextureFormatPixelBytes(img.pixelFormat);
        if (bytesPerPixel == 0u)
        {
            if (outError)
            {
                *outError = "invalid image pixel format";
            }
            return false;
        }
        if (img.mipLevels.empty())
        {
            if (outError)
            {
                *outError = "mip chain is empty";
            }
            return false;
        }
        for (size_t mipIndex = 0; mipIndex < img.mipLevels.size(); ++mipIndex)
        {
            const UdimMipImage &mip = img.mipLevels[mipIndex];
            const uint32_t expectedLevel = static_cast<uint32_t>(mipIndex);
            if (mip.level != expectedLevel)
            {
                if (outError)
                {
                    *outError = "mip levels must be contiguous starting at 0";
                }
                return false;
            }
            uint32_t expectedW = img.width;
            uint32_t expectedH = img.height;
            if (mipIndex > 0)
            {
                const UdimMipImage &prev = img.mipLevels[mipIndex - 1u];
                expectedW = std::max(1u, prev.width >> 1u);
                expectedH = std::max(1u, prev.height >> 1u);
            }
            const size_t expectedSamples =
                static_cast<size_t>(expectedW) * static_cast<size_t>(expectedH) *
                TextureFormatChannelCount(img.pixelFormat);
            if (mip.pixelFormat != img.pixelFormat || mip.width != expectedW || mip.height != expectedH ||
                mip.rgba.size() != expectedSamples)
            {
                if (outError)
                {
                    *outError = "provided mip level dimensions/pixels are invalid";
                }
                return false;
            }
        }
        const UdimMipImage &lastMip = img.mipLevels.back();
        if (lastMip.width != 1u || lastMip.height != 1u)
        {
            if (outError)
            {
                *outError = "mip chain must terminate at 1x1";
            }
            return false;
        }
    }

    return true;
}

bool ReadBytes(std::ifstream *stream, void *dst, size_t bytes)
{
    stream->read(reinterpret_cast<char *>(dst), static_cast<std::streamsize>(bytes));
    return stream->good();
}

} // namespace detail
} // namespace tilebin
} // namespace ybi
