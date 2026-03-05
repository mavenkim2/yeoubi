#include "tile_binary.h"
#include "tile_binary_detail.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <fstream>
#include <type_traits>

namespace ybi
{
namespace tilebin
{

static constexpr uint32_t kUdimMin = 1001;
static constexpr uint32_t kUdimMax = 1100;

static_assert(std::is_trivially_copyable<TileFileHeader>::value, "TileFileHeader must be POD");
static_assert(std::is_trivially_copyable<UdimEntry>::value, "UdimEntry must be POD");
static_assert(std::is_trivially_copyable<MipRecord>::value, "MipRecord must be POD");
static_assert(std::is_trivially_copyable<TileRecord>::value, "TileRecord must be POD");

bool WriteTileBinary(const std::filesystem::path &path,
                     int tileSize,
                     const std::vector<UdimImage> &images,
                     std::string *outError)
{
    if (!detail::ValidateImageBasics(images, tileSize, outError))
    {
        return false;
    }

    std::vector<UdimImage> sorted = images;
    std::sort(sorted.begin(), sorted.end(), [](const UdimImage &a, const UdimImage &b) { return a.udim < b.udim; });

    TileFileHeader header = {};
    std::memcpy(header.magic, "YBITILE3", 8);
    header.version = 3;
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
        std::vector<UdimMipImage> mipChain;
        if (!detail::BuildMipChain(img, &mipChain, outError))
        {
            return false;
        }

        UdimEntry &entry = entries[i];
        entry.udim = img.udim;
        entry.imageWidth = img.width;
        entry.imageHeight = img.height;
        entry.tileSize = static_cast<uint32_t>(tileSize);
        entry.mipCount = static_cast<uint32_t>(mipChain.size());
        entry.streamMipCount = 0u;
        entry.tailMipCount = 0u;
        for (const UdimMipImage &mip : mipChain)
        {
            if (std::max(mip.width, mip.height) >= static_cast<uint32_t>(tileSize))
            {
                ++entry.streamMipCount;
            }
            else
            {
                ++entry.tailMipCount;
            }
        }
        entry.mipRecordCount = entry.mipCount;

        entry.mipRecordOffset = static_cast<uint64_t>(out.tellp());
        std::vector<MipRecord> mipRecords(entry.mipRecordCount);
        out.write(reinterpret_cast<const char *>(mipRecords.data()),
                  static_cast<std::streamsize>(mipRecords.size() * sizeof(MipRecord)));
        if (!out.good())
        {
            if (outError)
            {
                *outError = "failed to write mip record placeholders: " + path.string();
            }
            return false;
        }

        entry.tileRecordOffset = static_cast<uint64_t>(out.tellp());
        entry.tileRecordCount = 0u;
        for (const UdimMipImage &mip : mipChain)
        {
            if (std::max(mip.width, mip.height) >= static_cast<uint32_t>(tileSize))
            {
                const uint32_t countX = static_cast<uint32_t>((mip.width + tileSize - 1) / tileSize);
                const uint32_t countY = static_cast<uint32_t>((mip.height + tileSize - 1) / tileSize);
                entry.tileRecordCount += countX * countY;
            }
        }

        std::vector<TileRecord> tileRecords(entry.tileRecordCount);
        out.write(reinterpret_cast<const char *>(tileRecords.data()),
                  static_cast<std::streamsize>(tileRecords.size() * sizeof(TileRecord)));
        if (!out.good())
        {
            if (outError)
            {
                *outError = "failed to write tile record placeholders: " + path.string();
            }
            return false;
        }

        entry.payloadOffset = static_cast<uint64_t>(out.tellp());
        uint32_t nextTileRecord = 0u;
        for (const UdimMipImage &mip : mipChain)
        {
            MipRecord mipRecord = {};
            mipRecord.mipLevel = mip.level;
            mipRecord.width = mip.width;
            mipRecord.height = mip.height;

            const bool isTail = (std::max(mip.width, mip.height) < static_cast<uint32_t>(tileSize));
            mipRecord.isTail = isTail ? 1u : 0u;
            if (isTail)
            {
                mipRecord.tileCountX = 1u;
                mipRecord.tileCountY = 1u;
                mipRecord.firstTileRecord = 0u;
                mipRecord.tileRecordCount = 0u;
                mipRecord.byteOffset = static_cast<uint64_t>(out.tellp());
                mipRecord.byteSize = static_cast<uint64_t>(mip.rgba.size() * sizeof(float));
                out.write(reinterpret_cast<const char *>(mip.rgba.data()),
                          static_cast<std::streamsize>(mipRecord.byteSize));
                if (!out.good())
                {
                    if (outError)
                    {
                        *outError = "failed writing tail mip payload: " + path.string();
                    }
                    return false;
                }
            }
            else
            {
                const uint32_t tileCountX = static_cast<uint32_t>((mip.width + tileSize - 1) / tileSize);
                const uint32_t tileCountY = static_cast<uint32_t>((mip.height + tileSize - 1) / tileSize);
                mipRecord.tileCountX = tileCountX;
                mipRecord.tileCountY = tileCountY;
                mipRecord.firstTileRecord = nextTileRecord;
                mipRecord.tileRecordCount = tileCountX * tileCountY;

                for (uint32_t ty = 0; ty < tileCountY; ++ty)
                {
                    for (uint32_t tx = 0; tx < tileCountX; ++tx)
                    {
                        const uint32_t tileIndex = nextTileRecord++;
                        int tileWidth = 0;
                        int tileHeight = 0;
                        detail::ExtractTileRgbaF32(mip.rgba,
                                                   static_cast<int>(mip.width),
                                                   static_cast<int>(mip.height),
                                                   static_cast<int>(tx),
                                                   static_cast<int>(ty),
                                                   tileSize,
                                                   tilePixels,
                                                   tileWidth,
                                                   tileHeight);

                        TileRecord &record = tileRecords[tileIndex];
                        record.mipLevel = mip.level;
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
            }

            mipRecords[mip.level] = mipRecord;
        }
        entry.payloadBytes = static_cast<uint64_t>(out.tellp()) - entry.payloadOffset;

        const std::streampos afterPayload = out.tellp();
        out.seekp(static_cast<std::streamoff>(entry.mipRecordOffset), std::ios::beg);
        out.write(reinterpret_cast<const char *>(mipRecords.data()),
                  static_cast<std::streamsize>(mipRecords.size() * sizeof(MipRecord)));
        out.seekp(static_cast<std::streamoff>(entry.tileRecordOffset), std::ios::beg);
        out.write(reinterpret_cast<const char *>(tileRecords.data()),
                  static_cast<std::streamsize>(tileRecords.size() * sizeof(TileRecord)));
        if (!out.good())
        {
            if (outError)
            {
                *outError = "failed rewriting records: " + path.string();
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

    if (outHeader.channels != 4 || outHeader.elementType != 1 || outHeader.udimCount == 0)
    {
        if (outError)
        {
            *outError = "invalid tile header data: " + path.string();
        }
        return false;
    }

    if (std::memcmp(outHeader.magic, "YBITILE2", 8) == 0 && outHeader.version == 2)
    {
        return detail::ReadTileBinaryV2(&in, path, outHeader, outEntries, outImages, outError);
    }

    if (std::memcmp(outHeader.magic, "YBITILE3", 8) != 0 || outHeader.version != 3)
    {
        if (outError)
        {
            *outError = "unsupported tile file version: " + path.string();
        }
        return false;
    }

    outEntries.assign(outHeader.udimCount, UdimEntry{});
    in.seekg(static_cast<std::streamoff>(outHeader.udimTableOffset), std::ios::beg);
    if (!detail::ReadBytes(&in, outEntries.data(), outEntries.size() * sizeof(UdimEntry)))
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
            entry.mipCount == 0 || entry.mipRecordCount != entry.mipCount)
        {
            if (outError)
            {
                *outError = "invalid UDIM entry in tile file: " + path.string();
            }
            return false;
        }

        std::vector<MipRecord> mipRecords(entry.mipRecordCount);
        in.seekg(static_cast<std::streamoff>(entry.mipRecordOffset), std::ios::beg);
        if (!detail::ReadBytes(&in, mipRecords.data(), mipRecords.size() * sizeof(MipRecord)))
        {
            if (outError)
            {
                *outError = "failed reading mip records: " + path.string();
            }
            return false;
        }

        std::vector<TileRecord> tileRecords(entry.tileRecordCount);
        in.seekg(static_cast<std::streamoff>(entry.tileRecordOffset), std::ios::beg);
        if (!detail::ReadBytes(&in, tileRecords.data(), tileRecords.size() * sizeof(TileRecord)))
        {
            if (outError)
            {
                *outError = "failed reading tile records: " + path.string();
            }
            return false;
        }

        const MipRecord &baseMip = mipRecords[0];
        if (baseMip.mipLevel != 0u || baseMip.width != entry.imageWidth || baseMip.height != entry.imageHeight)
        {
            if (outError)
            {
                *outError = "invalid base mip record in: " + path.string();
            }
            return false;
        }

        UdimImage image = {};
        image.udim = entry.udim;
        image.width = entry.imageWidth;
        image.height = entry.imageHeight;
        image.rgba.assign(static_cast<size_t>(entry.imageWidth) * static_cast<size_t>(entry.imageHeight) * 4u, 0.0f);

        if (baseMip.isTail != 0u)
        {
            const uint64_t expectedBytes =
                static_cast<uint64_t>(baseMip.width) * static_cast<uint64_t>(baseMip.height) * 4u * sizeof(float);
            if (baseMip.byteSize != expectedBytes)
            {
                if (outError)
                {
                    *outError = "base tail mip byte size mismatch in: " + path.string();
                }
                return false;
            }
            in.seekg(static_cast<std::streamoff>(baseMip.byteOffset), std::ios::beg);
            if (!detail::ReadBytes(&in, image.rgba.data(), static_cast<size_t>(baseMip.byteSize)))
            {
                if (outError)
                {
                    *outError = "failed reading base tail mip in: " + path.string();
                }
                return false;
            }
        }
        else
        {
            if (baseMip.firstTileRecord + baseMip.tileRecordCount > tileRecords.size())
            {
                if (outError)
                {
                    *outError = "base mip tile record range invalid in: " + path.string();
                }
                return false;
            }

            for (uint32_t i = 0u; i < baseMip.tileRecordCount; ++i)
            {
                const TileRecord &r = tileRecords[baseMip.firstTileRecord + i];
                if (r.mipLevel != 0u || r.width == 0 || r.height == 0 ||
                    r.width > entry.tileSize || r.height > entry.tileSize)
                {
                    if (outError)
                    {
                        *outError = "base tile record invalid in: " + path.string();
                    }
                    return false;
                }

                const uint64_t expectedBytes =
                    static_cast<uint64_t>(r.width) * static_cast<uint64_t>(r.height) * 4u * sizeof(float);
                if (r.byteSize != expectedBytes)
                {
                    if (outError)
                    {
                        *outError = "base tile byte size mismatch in: " + path.string();
                    }
                    return false;
                }

                std::vector<float> tile(static_cast<size_t>(r.width) * static_cast<size_t>(r.height) * 4u);
                in.seekg(static_cast<std::streamoff>(r.byteOffset), std::ios::beg);
                if (!detail::ReadBytes(&in, tile.data(), static_cast<size_t>(r.byteSize)))
                {
                    if (outError)
                    {
                        *outError = "failed reading base tile payload in: " + path.string();
                    }
                    return false;
                }

                const uint32_t x0 = r.tileX * entry.tileSize;
                const uint32_t y0 = r.tileY * entry.tileSize;
                if (x0 + r.width > entry.imageWidth || y0 + r.height > entry.imageHeight)
                {
                    if (outError)
                    {
                        *outError = "base tile bounds overflow in: " + path.string();
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
