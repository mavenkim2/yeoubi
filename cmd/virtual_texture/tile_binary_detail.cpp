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

void ExtractTileRgbaF32(const std::vector<float> &image,
                        int imageWidth,
                        int imageHeight,
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
    outTile.assign(static_cast<size_t>(outWidth) * static_cast<size_t>(outHeight) * 4u, 0.0f);
    for (int y = 0; y < outHeight; ++y)
    {
        const float *src = image.data() + (static_cast<size_t>(y0 + y) * static_cast<size_t>(imageWidth) +
                                           static_cast<size_t>(x0)) *
                                              4u;
        float *dst = outTile.data() + static_cast<size_t>(y) * static_cast<size_t>(outWidth) * 4u;
        std::memcpy(dst, src, static_cast<size_t>(outWidth) * 4u * sizeof(float));
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
            const size_t expectedPixels = static_cast<size_t>(expectedW) * static_cast<size_t>(expectedH) * 4u;
            if (mip.width != expectedW || mip.height != expectedH || mip.rgba.size() != expectedPixels)
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

bool ReadTileBinaryV2(std::ifstream *in,
                      const std::filesystem::path &path,
                      TileFileHeader &outHeader,
                      std::vector<UdimEntry> &outEntries,
                      std::vector<UdimImage> &outImages,
                      std::string *outError)
{
    std::vector<V2UdimEntry> v2Entries(outHeader.udimCount);
    in->seekg(static_cast<std::streamoff>(outHeader.udimTableOffset), std::ios::beg);
    if (!ReadBytes(in, v2Entries.data(), v2Entries.size() * sizeof(V2UdimEntry)))
    {
        if (outError)
        {
            *outError = "failed reading UDIM table: " + path.string();
        }
        return false;
    }

    outEntries.assign(v2Entries.size(), UdimEntry{});
    outImages.clear();
    outImages.reserve(v2Entries.size());

    for (size_t i = 0; i < v2Entries.size(); ++i)
    {
        const V2UdimEntry &entry = v2Entries[i];
        if (entry.udim < kUdimMin || entry.udim > kUdimMax || entry.tileSize == 0 ||
            entry.tileCount != entry.tileCountX * entry.tileCountY ||
            entry.tileRecordCount != entry.tileCount)
        {
            if (outError)
            {
                *outError = "invalid UDIM entry in tile file: " + path.string();
            }
            return false;
        }

        UdimEntry converted = {};
        converted.udim = entry.udim;
        converted.imageWidth = entry.imageWidth;
        converted.imageHeight = entry.imageHeight;
        converted.tileSize = entry.tileSize;
        converted.mipCount = 1u;
        converted.streamMipCount = 1u;
        converted.tailMipCount = 0u;
        converted.tileRecordOffset = entry.tileRecordOffset;
        converted.tileRecordCount = entry.tileRecordCount;
        converted.payloadOffset = entry.payloadOffset;
        converted.payloadBytes = entry.payloadBytes;
        outEntries[i] = converted;

        std::vector<V2TileRecord> records(entry.tileRecordCount);
        in->seekg(static_cast<std::streamoff>(entry.tileRecordOffset), std::ios::beg);
        if (!ReadBytes(in, records.data(), records.size() * sizeof(V2TileRecord)))
        {
            if (outError)
            {
                *outError = "failed reading tile records: " + path.string();
            }
            return false;
        }

        UdimImage image = {};
        image.udim = entry.udim;
        image.width = entry.imageWidth;
        image.height = entry.imageHeight;
        image.mipLevels.resize(1u);
        image.mipLevels[0].level = 0u;
        image.mipLevels[0].width = entry.imageWidth;
        image.mipLevels[0].height = entry.imageHeight;
        image.mipLevels[0].rgba.assign(static_cast<size_t>(entry.imageWidth) * static_cast<size_t>(entry.imageHeight) * 4u,
                                       0.0f);

        for (const V2TileRecord &r : records)
        {
            if (r.tileX >= entry.tileCountX || r.tileY >= entry.tileCountY ||
                r.width == 0 || r.height == 0 || r.width > entry.tileSize || r.height > entry.tileSize)
            {
                if (outError)
                {
                    *outError = "tile record out of range in: " + path.string();
                }
                return false;
            }

            const uint64_t expectedBytes =
                static_cast<uint64_t>(r.width) * static_cast<uint64_t>(r.height) * 4u * sizeof(float);
            if (r.byteSize != expectedBytes)
            {
                if (outError)
                {
                    *outError = "tile byte size mismatch in: " + path.string();
                }
                return false;
            }

            std::vector<float> tile(static_cast<size_t>(r.width) * static_cast<size_t>(r.height) * 4u);
            in->seekg(static_cast<std::streamoff>(r.byteOffset), std::ios::beg);
            if (!ReadBytes(in, tile.data(), static_cast<size_t>(r.byteSize)))
            {
                if (outError)
                {
                    *outError = "failed reading tile payload in: " + path.string();
                }
                return false;
            }

            const uint32_t x0 = r.tileX * entry.tileSize;
            const uint32_t y0 = r.tileY * entry.tileSize;
            if (x0 + r.width > entry.imageWidth || y0 + r.height > entry.imageHeight)
            {
                if (outError)
                {
                    *outError = "tile bounds overflow in: " + path.string();
                }
                return false;
            }

            for (uint32_t y = 0; y < r.height; ++y)
            {
                const float *src = tile.data() + static_cast<size_t>(y) * static_cast<size_t>(r.width) * 4u;
                float *dst = image.mipLevels[0].rgba.data() +
                             (static_cast<size_t>(y0 + y) * static_cast<size_t>(entry.imageWidth) + static_cast<size_t>(x0)) *
                                 4u;
                std::memcpy(dst, src, static_cast<size_t>(r.width) * 4u * sizeof(float));
            }
        }

        outImages.push_back(std::move(image));
    }

    return true;
}

} // namespace detail
} // namespace tilebin
} // namespace ybi
