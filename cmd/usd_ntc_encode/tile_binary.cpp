#include "tile_binary.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <fstream>
#include <type_traits>

namespace ybi
{
namespace tilebin
{

static_assert(std::is_trivially_copyable<TileFileHeader>::value, "TileFileHeader must be POD");
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

} // namespace

bool WriteTileBinary(const std::filesystem::path &path,
                     int imageWidth,
                     int imageHeight,
                     int tileSize,
                     const std::vector<float> &image,
                     std::string *outError)
{
    const uint32_t tilesX = static_cast<uint32_t>((imageWidth + tileSize - 1) / tileSize);
    const uint32_t tilesY = static_cast<uint32_t>((imageHeight + tileSize - 1) / tileSize);
    const uint32_t tileCount = tilesX * tilesY;

    TileFileHeader header = {};
    std::memcpy(header.magic, "YBITILE1", 8);
    header.tileSize = static_cast<uint32_t>(tileSize);
    header.imageWidth = static_cast<uint32_t>(imageWidth);
    header.imageHeight = static_cast<uint32_t>(imageHeight);
    header.tileCountX = tilesX;
    header.tileCountY = tilesY;
    header.tileCount = tileCount;

    std::vector<TileRecord> records(tileCount);
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
    out.write(reinterpret_cast<const char *>(records.data()),
              static_cast<std::streamsize>(records.size() * sizeof(TileRecord)));
    if (!out.good())
    {
        if (outError)
        {
            *outError = "failed to write tile file header: " + path.string();
        }
        return false;
    }

    std::vector<float> tilePixels;
    for (uint32_t ty = 0; ty < tilesY; ++ty)
    {
        for (uint32_t tx = 0; tx < tilesX; ++tx)
        {
            const uint32_t tileIndex = ty * tilesX + tx;
            int tileWidth = 0;
            int tileHeight = 0;
            ExtractTileRgbaF32(
                image, imageWidth, imageHeight, static_cast<int>(tx), static_cast<int>(ty), tileSize, tilePixels, tileWidth, tileHeight);

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

    out.seekp(0, std::ios::beg);
    out.write(reinterpret_cast<const char *>(&header), sizeof(header));
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
    out.close();
    return true;
}

bool ReadTileBinary(const std::filesystem::path &path,
                    TileFileHeader &outHeader,
                    std::vector<TileRecord> &outRecords,
                    std::vector<float> &outImage,
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
    if (std::memcmp(outHeader.magic, "YBITILE1", 8) != 0)
    {
        if (outError)
        {
            *outError = "bad tile magic: " + path.string();
        }
        return false;
    }
    if (outHeader.channels != 4 || outHeader.elementType != 1 || outHeader.tileSize == 0)
    {
        if (outError)
        {
            *outError = "unsupported tile format in: " + path.string();
        }
        return false;
    }
    if (outHeader.tileCount != outHeader.tileCountX * outHeader.tileCountY)
    {
        if (outError)
        {
            *outError = "tile count mismatch in: " + path.string();
        }
        return false;
    }

    outRecords.resize(outHeader.tileCount);
    in.read(reinterpret_cast<char *>(outRecords.data()),
            static_cast<std::streamsize>(outRecords.size() * sizeof(TileRecord)));
    if (!in.good())
    {
        if (outError)
        {
            *outError = "failed reading tile records: " + path.string();
        }
        return false;
    }

    const size_t width = static_cast<size_t>(outHeader.imageWidth);
    const size_t height = static_cast<size_t>(outHeader.imageHeight);
    outImage.assign(width * height * 4u, 0.0f);

    for (const TileRecord &r : outRecords)
    {
        if (r.tileX >= outHeader.tileCountX || r.tileY >= outHeader.tileCountY)
        {
            if (outError)
            {
                *outError = "tile record out of range in: " + path.string();
            }
            return false;
        }
        if (r.width == 0 || r.height == 0 || r.width > outHeader.tileSize || r.height > outHeader.tileSize)
        {
            if (outError)
            {
                *outError = "tile record size invalid in: " + path.string();
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

        const uint32_t x0 = r.tileX * outHeader.tileSize;
        const uint32_t y0 = r.tileY * outHeader.tileSize;
        if (x0 + r.width > outHeader.imageWidth || y0 + r.height > outHeader.imageHeight)
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
            float *dst = outImage.data() +
                         (static_cast<size_t>(y0 + y) * static_cast<size_t>(outHeader.imageWidth) + static_cast<size_t>(x0)) * 4u;
            std::memcpy(dst, src, static_cast<size_t>(r.width) * 4u * sizeof(float));
        }
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
