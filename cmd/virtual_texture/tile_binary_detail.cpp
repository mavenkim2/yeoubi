#include "tile_binary_detail.h"

#include <algorithm>
#include <cstring>
#include <unordered_set>

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
                (static_cast<size_t>(y) * static_cast<size_t>(*dstW) + static_cast<size_t>(x)) * 4u;
            for (uint32_t c = 0; c < 4u; ++c)
            {
                const size_t i00 =
                    (static_cast<size_t>(sy0) * static_cast<size_t>(srcW) + static_cast<size_t>(sx0)) * 4u + c;
                const size_t i10 =
                    (static_cast<size_t>(sy0) * static_cast<size_t>(srcW) + static_cast<size_t>(sx1)) * 4u + c;
                const size_t i01 =
                    (static_cast<size_t>(sy1) * static_cast<size_t>(srcW) + static_cast<size_t>(sx0)) * 4u + c;
                const size_t i11 =
                    (static_cast<size_t>(sy1) * static_cast<size_t>(srcW) + static_cast<size_t>(sx1)) * 4u + c;
                (*dst)[outBase + c] =
                    0.25f * (src[i00] + src[i10] + src[i01] + src[i11]);
            }
        }
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

    std::unordered_set<uint32_t> seen;
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
        if (!seen.insert(img.udim).second)
        {
            if (outError)
            {
                *outError = "duplicate UDIM id";
            }
            return false;
        }
        if (img.width == 0 || img.height == 0)
        {
            if (outError)
            {
                *outError = "invalid image dimensions";
            }
            return false;
        }
        const size_t expected = static_cast<size_t>(img.width) * static_cast<size_t>(img.height) * 4u;
        if (img.rgba.size() != expected)
        {
            if (outError)
            {
                *outError = "image pixel count mismatch";
            }
            return false;
        }
    }
    return true;
}

bool BuildMipChain(const UdimImage &img,
                   std::vector<UdimMipImage> *outMipChain,
                   std::string *outError)
{
    outMipChain->clear();

    UdimMipImage mip0 = {};
    mip0.level = 0u;
    mip0.width = img.width;
    mip0.height = img.height;
    mip0.rgba = img.rgba;
    outMipChain->push_back(std::move(mip0));

    if (!img.mipLevels.empty())
    {
        std::vector<UdimMipImage> provided = img.mipLevels;
        std::sort(provided.begin(), provided.end(),
                  [](const UdimMipImage &a, const UdimMipImage &b) { return a.level < b.level; });

        size_t next = 0;
        if (!provided.empty() && provided[0].level == 0u)
        {
            const UdimMipImage &src0 = provided[0];
            const size_t expected = static_cast<size_t>(img.width) * static_cast<size_t>(img.height) * 4u;
            if (src0.width != img.width || src0.height != img.height || src0.rgba.size() != expected)
            {
                if (outError)
                {
                    *outError = "mip level 0 does not match base image";
                }
                return false;
            }
            next = 1;
        }

        uint32_t expectedLevel = 1u;
        while (next < provided.size())
        {
            const UdimMipImage &src = provided[next];
            if (src.level != expectedLevel)
            {
                if (outError)
                {
                    *outError = "mip levels must be contiguous starting at 0";
                }
                return false;
            }

            const uint32_t expectedW =
                std::max(1u, outMipChain->back().width >> 1u);
            const uint32_t expectedH =
                std::max(1u, outMipChain->back().height >> 1u);
            const size_t expectedPixels = static_cast<size_t>(src.width) * static_cast<size_t>(src.height) * 4u;
            if (src.width != expectedW || src.height != expectedH || src.rgba.size() != expectedPixels)
            {
                if (outError)
                {
                    *outError = "provided mip level dimensions/pixels are invalid";
                }
                return false;
            }

            outMipChain->push_back(src);
            ++next;
            ++expectedLevel;
        }
    }

    while (outMipChain->back().width > 1u || outMipChain->back().height > 1u)
    {
        UdimMipImage next = {};
        next.level = static_cast<uint32_t>(outMipChain->size());
        Downsample2x2(outMipChain->back().rgba,
                      outMipChain->back().width,
                      outMipChain->back().height,
                      &next.rgba,
                      &next.width,
                      &next.height);
        outMipChain->push_back(std::move(next));
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
        image.rgba.assign(static_cast<size_t>(entry.imageWidth) * static_cast<size_t>(entry.imageHeight) * 4u, 0.0f);

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
                float *dst = image.rgba.data() +
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
