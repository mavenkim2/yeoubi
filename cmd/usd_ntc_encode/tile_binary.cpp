#include "tile_binary.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <fstream>
#include <type_traits>
#include <unordered_set>

namespace ybi
{
namespace tilebin
{

static constexpr uint32_t kUdimMin = 1001;
static constexpr uint32_t kUdimMax = 1100;

static_assert(std::is_trivially_copyable<TileFileHeader>::value, "TileFileHeader must be POD");
static_assert(std::is_trivially_copyable<UdimEntry>::value, "UdimEntry must be POD");
static_assert(std::is_trivially_copyable<TileRecord>::value, "TileRecord must be POD");

namespace
{

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

bool ValidateImages(const std::vector<UdimImage> &images, int tileSize, std::string *outError)
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

} // namespace

bool WriteTileBinary(const std::filesystem::path &path,
                     int tileSize,
                     const std::vector<UdimImage> &images,
                     std::string *outError)
{
    if (!ValidateImages(images, tileSize, outError))
    {
        return false;
    }

    std::vector<UdimImage> sorted = images;
    std::sort(sorted.begin(), sorted.end(), [](const UdimImage &a, const UdimImage &b) { return a.udim < b.udim; });

    TileFileHeader header = {};
    std::memcpy(header.magic, "YBITILE2", 8);
    header.version = 2;
    header.channels = 4;
    header.elementType = 1;
    header.udimCount = static_cast<uint32_t>(sorted.size());
    header.udimTableOffset = sizeof(TileFileHeader);

    std::vector<UdimEntry> entries(sorted.size());

    std::fstream out(path, std::ios::binary | std::ios::out | std::ios::trunc);
    if (!out.is_open())
    {
        if (outError)
        {
            *outError = "failed to open output tile file: " + path.string();
        }
        return false;
    }

    out.write(reinterpret_cast<const char *>(&header), sizeof(header));
    out.write(reinterpret_cast<const char *>(entries.data()),
              static_cast<std::streamsize>(entries.size() * sizeof(UdimEntry)));
    if (!out.good())
    {
        if (outError)
        {
            *outError = "failed to write header/table placeholders: " + path.string();
        }
        return false;
    }

    std::vector<float> tilePixels;
    for (size_t i = 0; i < sorted.size(); ++i)
    {
        const UdimImage &img = sorted[i];
        UdimEntry &entry = entries[i];
        entry.udim = img.udim;
        entry.imageWidth = img.width;
        entry.imageHeight = img.height;
        entry.tileSize = static_cast<uint32_t>(tileSize);
        entry.tileCountX = static_cast<uint32_t>((img.width + tileSize - 1) / tileSize);
        entry.tileCountY = static_cast<uint32_t>((img.height + tileSize - 1) / tileSize);
        entry.tileCount = entry.tileCountX * entry.tileCountY;
        entry.tileRecordCount = entry.tileCount;

        entry.tileRecordOffset = static_cast<uint64_t>(out.tellp());
        std::vector<TileRecord> records(entry.tileCount);
        out.write(reinterpret_cast<const char *>(records.data()),
                  static_cast<std::streamsize>(records.size() * sizeof(TileRecord)));
        if (!out.good())
        {
            if (outError)
            {
                *outError = "failed to write tile record placeholders: " + path.string();
            }
            return false;
        }

        entry.payloadOffset = static_cast<uint64_t>(out.tellp());
        for (uint32_t ty = 0; ty < entry.tileCountY; ++ty)
        {
            for (uint32_t tx = 0; tx < entry.tileCountX; ++tx)
            {
                const uint32_t tileIndex = ty * entry.tileCountX + tx;
                int tileWidth = 0;
                int tileHeight = 0;
                ExtractTileRgbaF32(img.rgba,
                                   static_cast<int>(img.width),
                                   static_cast<int>(img.height),
                                   static_cast<int>(tx),
                                   static_cast<int>(ty),
                                   tileSize,
                                   tilePixels,
                                   tileWidth,
                                   tileHeight);

                TileRecord &record = records[tileIndex];
                record.tileX = tx;
                record.tileY = ty;
                record.width = static_cast<uint32_t>(tileWidth);
                record.height = static_cast<uint32_t>(tileHeight);
                record.byteOffset = static_cast<uint64_t>(out.tellp());
                record.byteSize = static_cast<uint64_t>(tilePixels.size() * sizeof(float));

                out.write(reinterpret_cast<const char *>(tilePixels.data()),
                          static_cast<std::streamsize>(record.byteSize));
                if (!out.good())
                {
                    if (outError)
                    {
                        *outError = "failed writing tile payload: " + path.string();
                    }
                    return false;
                }
            }
        }
        entry.payloadBytes = static_cast<uint64_t>(out.tellp()) - entry.payloadOffset;

        const std::streampos afterPayload = out.tellp();
        out.seekp(static_cast<std::streamoff>(entry.tileRecordOffset), std::ios::beg);
        out.write(reinterpret_cast<const char *>(records.data()),
                  static_cast<std::streamsize>(records.size() * sizeof(TileRecord)));
        if (!out.good())
        {
            if (outError)
            {
                *outError = "failed rewriting tile records: " + path.string();
            }
            return false;
        }
        out.seekp(afterPayload, std::ios::beg);
    }

    out.seekp(0, std::ios::beg);
    out.write(reinterpret_cast<const char *>(&header), sizeof(header));
    out.write(reinterpret_cast<const char *>(entries.data()),
              static_cast<std::streamsize>(entries.size() * sizeof(UdimEntry)));
    if (!out.good())
    {
        if (outError)
        {
            *outError = "failed rewriting header/table: " + path.string();
        }
        return false;
    }

    out.close();
    return true;
}

bool ReadTileBinary(const std::filesystem::path &path,
                    TileFileHeader &outHeader,
                    std::vector<UdimEntry> &outEntries,
                    std::vector<UdimImage> &outImages,
                    std::string *outError)
{
    std::ifstream in(path, std::ios::binary);
    if (!in.is_open())
    {
        if (outError)
        {
            *outError = "failed to open tile file: " + path.string();
        }
        return false;
    }

    in.read(reinterpret_cast<char *>(&outHeader), sizeof(outHeader));
    if (!in.good())
    {
        if (outError)
        {
            *outError = "failed reading tile header: " + path.string();
        }
        return false;
    }
    if (std::memcmp(outHeader.magic, "YBITILE2", 8) != 0 || outHeader.version != 2)
    {
        if (outError)
        {
            *outError = "unsupported tile file version: " + path.string();
        }
        return false;
    }
    if (outHeader.channels != 4 || outHeader.elementType != 1 || outHeader.udimCount == 0)
    {
        if (outError)
        {
            *outError = "invalid tile header data: " + path.string();
        }
        return false;
    }

    outEntries.assign(outHeader.udimCount, UdimEntry{});
    in.seekg(static_cast<std::streamoff>(outHeader.udimTableOffset), std::ios::beg);
    in.read(reinterpret_cast<char *>(outEntries.data()),
            static_cast<std::streamsize>(outEntries.size() * sizeof(UdimEntry)));
    if (!in.good())
    {
        if (outError)
        {
            *outError = "failed reading UDIM table: " + path.string();
        }
        return false;
    }

    outImages.clear();
    outImages.reserve(outEntries.size());
    for (const UdimEntry &entry : outEntries)
    {
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

        std::vector<TileRecord> records(entry.tileRecordCount);
        in.seekg(static_cast<std::streamoff>(entry.tileRecordOffset), std::ios::beg);
        in.read(reinterpret_cast<char *>(records.data()),
                static_cast<std::streamsize>(records.size() * sizeof(TileRecord)));
        if (!in.good())
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

        for (const TileRecord &r : records)
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

            const uint64_t expectedBytes = static_cast<uint64_t>(r.width) * static_cast<uint64_t>(r.height) * 4u * sizeof(float);
            if (r.byteSize != expectedBytes)
            {
                if (outError)
                {
                    *outError = "tile byte size mismatch in: " + path.string();
                }
                return false;
            }

            std::vector<float> tile(static_cast<size_t>(r.width) * static_cast<size_t>(r.height) * 4u);
            in.seekg(static_cast<std::streamoff>(r.byteOffset), std::ios::beg);
            in.read(reinterpret_cast<char *>(tile.data()), static_cast<std::streamsize>(r.byteSize));
            if (!in.good())
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
                             (static_cast<size_t>(y0 + y) * static_cast<size_t>(entry.imageWidth) + static_cast<size_t>(x0)) * 4u;
                std::memcpy(dst, src, static_cast<size_t>(r.width) * 4u * sizeof(float));
            }
        }

        outImages.push_back(std::move(image));
    }

    return true;
}

bool DiffImagesExact(const std::vector<float> &a, const std::vector<float> &b, float eps, DiffStats *outStats)
{
    DiffStats stats = {};
    stats.firstMismatch = a.size();
    if (a.size() != b.size())
    {
        stats.mismatchCount = 1;
        stats.firstMismatch = 0;
        if (outStats)
        {
            *outStats = stats;
        }
        return false;
    }

    double sumSq = 0.0;
    double sumAbs = 0.0;
    for (size_t i = 0; i < a.size(); ++i)
    {
        const float diff = std::fabs(a[i] - b[i]);
        stats.maxAbs = std::max(stats.maxAbs, diff);
        sumAbs += diff;
        sumSq += static_cast<double>(diff) * static_cast<double>(diff);
        if (diff > eps)
        {
            ++stats.mismatchCount;
            if (stats.firstMismatch == a.size())
            {
                stats.firstMismatch = i;
            }
        }
    }

    const double denom = a.empty() ? 1.0 : static_cast<double>(a.size());
    stats.meanAbs = sumAbs / denom;
    stats.rmse = std::sqrt(sumSq / denom);
    if (outStats)
    {
        *outStats = stats;
    }
    return stats.mismatchCount == 0;
}

} // namespace tilebin
} // namespace ybi
