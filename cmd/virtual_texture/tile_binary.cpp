#include "tile_binary.h"
#include "tile_binary_detail.h"

#include "miniz.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <fstream>
#include <limits>
#include <type_traits>

namespace ybi
{
namespace tilebin
{

static constexpr uint32_t kUdimMin = 1001;
static constexpr uint32_t kUdimMax = 1100;
static constexpr uint32_t kCompressionNone = 0u;
static constexpr uint32_t kCompressionDeflate = 1u;

static_assert(std::is_trivially_copyable<TileFileHeader>::value, "TileFileHeader must be POD");
static_assert(std::is_trivially_copyable<UdimEntry>::value, "UdimEntry must be POD");
static_assert(std::is_trivially_copyable<MipRecord>::value, "MipRecord must be POD");
static_assert(std::is_trivially_copyable<TileRecord>::value, "TileRecord must be POD");

namespace
{

struct V4MipRecord
{
    uint32_t mipLevel = 0;
    uint32_t width = 0;
    uint32_t height = 0;
    uint32_t tileCountX = 0;
    uint32_t tileCountY = 0;
    uint32_t firstTileRecord = 0;
    uint32_t tileRecordCount = 0;
    uint32_t isTail = 0;
    uint64_t byteOffset = 0;
    uint64_t storedByteSize = 0;
    uint64_t rawByteSize = 0;
    uint32_t compression = 0;
    uint32_t reserved0 = 0;
};

struct V4TileRecord
{
    uint32_t mipLevel = 0;
    uint32_t tileX = 0;
    uint32_t tileY = 0;
    uint32_t width = 0;
    uint32_t height = 0;
    uint32_t compression = 0;
    uint32_t reserved0 = 0;
    uint64_t byteOffset = 0;
    uint64_t storedByteSize = 0;
    uint64_t rawByteSize = 0;
};

struct ReadMipRecord
{
    uint32_t mipLevel = 0;
    uint32_t width = 0;
    uint32_t height = 0;
    uint32_t tileCountX = 0;
    uint32_t tileCountY = 0;
    uint32_t firstTileRecord = 0;
    uint32_t tileRecordCount = 0;
    uint32_t isTail = 0;
    uint64_t byteOffset = 0;
    uint64_t storedByteSize = 0;
    uint64_t rawByteSize = 0;
    uint32_t compression = 0;
};

struct ReadTileRecord
{
    uint32_t mipLevel = 0;
    uint32_t tileX = 0;
    uint32_t tileY = 0;
    uint32_t width = 0;
    uint32_t height = 0;
    uint64_t byteOffset = 0;
    uint64_t storedByteSize = 0;
    uint64_t rawByteSize = 0;
    uint32_t compression = 0;
};

static_assert(std::is_trivially_copyable<V4MipRecord>::value, "V4MipRecord must be POD");
static_assert(std::is_trivially_copyable<V4TileRecord>::value, "V4TileRecord must be POD");

bool EnsureMzUlongRange(uint64_t bytes, const std::filesystem::path &path, std::string *outError)
{
    if (bytes > static_cast<uint64_t>(std::numeric_limits<mz_ulong>::max()))
    {
        if (outError)
        {
            *outError = "payload too large for miniz API: " + path.string();
        }
        return false;
    }
    return true;
}

bool CompressPayload(const unsigned char *rawBytes,
                     uint64_t rawByteCount,
                     std::vector<unsigned char> *outStored,
                     uint32_t *outCompression,
                     const std::filesystem::path &path,
                     std::string *outError)
{
    if (!EnsureMzUlongRange(rawByteCount, path, outError))
    {
        return false;
    }

    if (rawByteCount == 0u)
    {
        outStored->clear();
        *outCompression = kCompressionNone;
        return true;
    }

    const mz_ulong bound = mz_compressBound(static_cast<mz_ulong>(rawByteCount));
    std::vector<unsigned char> compressed(static_cast<size_t>(bound));
    mz_ulong compressedBytes = bound;
    const int rc = mz_compress2(compressed.data(),
                                &compressedBytes,
                                rawBytes,
                                static_cast<mz_ulong>(rawByteCount),
                                MZ_DEFAULT_LEVEL);
    if (rc == MZ_OK && compressedBytes < static_cast<mz_ulong>(rawByteCount))
    {
        compressed.resize(static_cast<size_t>(compressedBytes));
        *outStored = std::move(compressed);
        *outCompression = kCompressionDeflate;
        return true;
    }

    outStored->assign(rawBytes, rawBytes + rawByteCount);
    *outCompression = kCompressionNone;
    return true;
}

bool DecodePayload(const std::vector<unsigned char> &stored,
                   uint64_t rawByteCount,
                   uint32_t compression,
                   std::vector<unsigned char> *outRaw,
                   const std::filesystem::path &path,
                   std::string *outError)
{
    if (!EnsureMzUlongRange(rawByteCount, path, outError))
    {
        return false;
    }
    if (!EnsureMzUlongRange(stored.size(), path, outError))
    {
        return false;
    }

    if (compression == kCompressionNone)
    {
        if (stored.size() != static_cast<size_t>(rawByteCount))
        {
            if (outError)
            {
                *outError = "raw payload size mismatch in: " + path.string();
            }
            return false;
        }
        *outRaw = stored;
        return true;
    }

    if (compression != kCompressionDeflate)
    {
        if (outError)
        {
            *outError = "unsupported tile payload compression in: " + path.string();
        }
        return false;
    }

    outRaw->assign(static_cast<size_t>(rawByteCount), 0u);
    mz_ulong decodedBytes = static_cast<mz_ulong>(rawByteCount);
    const int rc = mz_uncompress(outRaw->data(),
                                 &decodedBytes,
                                 stored.data(),
                                 static_cast<mz_ulong>(stored.size()));
    if (rc != MZ_OK || decodedBytes != static_cast<mz_ulong>(rawByteCount))
    {
        if (outError)
        {
            *outError = "failed to decompress tile payload in: " + path.string();
        }
        return false;
    }
    return true;
}

bool ReadDecodedPayload(std::ifstream *in,
                        uint64_t byteOffset,
                        uint64_t storedByteSize,
                        uint64_t rawByteSize,
                        uint32_t compression,
                        std::vector<unsigned char> *outRaw,
                        const std::filesystem::path &path,
                        std::string *outError)
{
    std::vector<unsigned char> stored(static_cast<size_t>(storedByteSize));
    in->seekg(static_cast<std::streamoff>(byteOffset), std::ios::beg);
    if (!detail::ReadBytes(in, stored.data(), stored.size()))
    {
        if (outError)
        {
            *outError = "failed reading compressed payload in: " + path.string();
        }
        return false;
    }
    return DecodePayload(stored, rawByteSize, compression, outRaw, path, outError);
}

bool ReconstructBaseImageFromRecords(const std::filesystem::path &path,
                                     const UdimEntry &entry,
                                     const std::vector<ReadMipRecord> &mipRecords,
                                     const std::vector<ReadTileRecord> &tileRecords,
                                     std::ifstream *in,
                                     UdimImage *outImage,
                                     std::string *outError)
{
    const ReadMipRecord &baseMip = mipRecords[0];
    if (baseMip.mipLevel != 0u || baseMip.width != entry.imageWidth || baseMip.height != entry.imageHeight)
    {
        if (outError)
        {
            *outError = "invalid base mip record in: " + path.string();
        }
        return false;
    }

    outImage->udim = entry.udim;
    outImage->width = entry.imageWidth;
    outImage->height = entry.imageHeight;
    outImage->rgba.assign(static_cast<size_t>(entry.imageWidth) * static_cast<size_t>(entry.imageHeight) * 4u, 0.0f);

    if (baseMip.isTail != 0u)
    {
        const uint64_t expectedBytes =
            static_cast<uint64_t>(baseMip.width) * static_cast<uint64_t>(baseMip.height) * 4u * sizeof(float);
        if (baseMip.rawByteSize != expectedBytes)
        {
            if (outError)
            {
                *outError = "base tail mip byte size mismatch in: " + path.string();
            }
            return false;
        }

        std::vector<unsigned char> rawPayload;
        if (!ReadDecodedPayload(in,
                                baseMip.byteOffset,
                                baseMip.storedByteSize,
                                baseMip.rawByteSize,
                                baseMip.compression,
                                &rawPayload,
                                path,
                                outError))
        {
            return false;
        }

        std::memcpy(outImage->rgba.data(), rawPayload.data(), static_cast<size_t>(baseMip.rawByteSize));
        return true;
    }

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
        const ReadTileRecord &r = tileRecords[baseMip.firstTileRecord + i];
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
        if (r.rawByteSize != expectedBytes)
        {
            if (outError)
            {
                *outError = "base tile byte size mismatch in: " + path.string();
            }
            return false;
        }

        std::vector<unsigned char> rawPayload;
        if (!ReadDecodedPayload(in,
                                r.byteOffset,
                                r.storedByteSize,
                                r.rawByteSize,
                                r.compression,
                                &rawPayload,
                                path,
                                outError))
        {
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

        const float *tileF32 = reinterpret_cast<const float *>(rawPayload.data());
        for (uint32_t y = 0; y < r.height; ++y)
        {
            const float *src = tileF32 + static_cast<size_t>(y) * static_cast<size_t>(r.width) * 4u;
            float *dst = outImage->rgba.data() +
                         (static_cast<size_t>(y0 + y) * static_cast<size_t>(entry.imageWidth) + static_cast<size_t>(x0)) *
                             4u;
            std::memcpy(dst, src, static_cast<size_t>(r.width) * 4u * sizeof(float));
        }
    }

    return true;
}

bool ReadEntries(std::ifstream *in,
                 const TileFileHeader &header,
                 std::vector<UdimEntry> *outEntries,
                 const std::filesystem::path &path,
                 std::string *outError)
{
    outEntries->assign(header.udimCount, UdimEntry{});
    in->seekg(static_cast<std::streamoff>(header.udimTableOffset), std::ios::beg);
    if (!detail::ReadBytes(in, outEntries->data(), outEntries->size() * sizeof(UdimEntry)))
    {
        if (outError)
        {
            *outError = "failed reading UDIM table: " + path.string();
        }
        return false;
    }
    return true;
}

} // namespace

bool WriteTileBinary(const std::filesystem::path &path,
                     int tileSize,
                     const std::vector<UdimImage> &images,
                     std::string *outError)
{
    if (!detail::ValidateImageBasics(images, tileSize, outError))
    {
        return false;
    }

    TileFileHeader header = {};
    std::memcpy(header.magic, "YBITILE4", 8);
    header.version = 4;
    header.channels = 4;
    header.elementType = 1;
    header.udimCount = static_cast<uint32_t>(images.size());
    header.udimTableOffset = sizeof(TileFileHeader);

    std::vector<UdimEntry> entries(images.size());

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
    for (size_t i = 0; i < images.size(); ++i)
    {
        const UdimImage &img = images[i];
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
        std::vector<V4MipRecord> mipRecords(entry.mipRecordCount);
        out.write(reinterpret_cast<const char *>(mipRecords.data()),
                  static_cast<std::streamsize>(mipRecords.size() * sizeof(V4MipRecord)));
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

        std::vector<V4TileRecord> tileRecords(entry.tileRecordCount);
        out.write(reinterpret_cast<const char *>(tileRecords.data()),
                  static_cast<std::streamsize>(tileRecords.size() * sizeof(V4TileRecord)));
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
            V4MipRecord mipRecord = {};
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
                mipRecord.rawByteSize = static_cast<uint64_t>(mip.rgba.size() * sizeof(float));

                std::vector<unsigned char> compressed;
                if (!CompressPayload(reinterpret_cast<const unsigned char *>(mip.rgba.data()),
                                     mipRecord.rawByteSize,
                                     &compressed,
                                     &mipRecord.compression,
                                     path,
                                     outError))
                {
                    return false;
                }
                mipRecord.storedByteSize = static_cast<uint64_t>(compressed.size());
                out.write(reinterpret_cast<const char *>(compressed.data()),
                          static_cast<std::streamsize>(mipRecord.storedByteSize));
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

                        V4TileRecord &record = tileRecords[tileIndex];
                        record.mipLevel = mip.level;
                        record.tileX = tx;
                        record.tileY = ty;
                        record.width = static_cast<uint32_t>(tileWidth);
                        record.height = static_cast<uint32_t>(tileHeight);
                        record.byteOffset = static_cast<uint64_t>(out.tellp());
                        record.rawByteSize = static_cast<uint64_t>(tilePixels.size() * sizeof(float));

                        std::vector<unsigned char> compressed;
                        if (!CompressPayload(reinterpret_cast<const unsigned char *>(tilePixels.data()),
                                             record.rawByteSize,
                                             &compressed,
                                             &record.compression,
                                             path,
                                             outError))
                        {
                            return false;
                        }
                        record.storedByteSize = static_cast<uint64_t>(compressed.size());

                        out.write(reinterpret_cast<const char *>(compressed.data()),
                                  static_cast<std::streamsize>(record.storedByteSize));
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
                  static_cast<std::streamsize>(mipRecords.size() * sizeof(V4MipRecord)));
        out.seekp(static_cast<std::streamoff>(entry.tileRecordOffset), std::ios::beg);
        out.write(reinterpret_cast<const char *>(tileRecords.data()),
                  static_cast<std::streamsize>(tileRecords.size() * sizeof(V4TileRecord)));
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

    const bool isV3 = (std::memcmp(outHeader.magic, "YBITILE3", 8) == 0 && outHeader.version == 3);
    const bool isV4 = (std::memcmp(outHeader.magic, "YBITILE4", 8) == 0 && outHeader.version == 4);
    if (!isV3 && !isV4)
    {
        if (outError)
        {
            *outError = "unsupported tile file version: " + path.string();
        }
        return false;
    }

    if (!ReadEntries(&in, outHeader, &outEntries, path, outError))
    {
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

        std::vector<ReadMipRecord> readMips;
        readMips.resize(entry.mipRecordCount);

        std::vector<ReadTileRecord> readTiles;
        readTiles.resize(entry.tileRecordCount);

        if (isV3)
        {
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
            for (size_t i = 0; i < mipRecords.size(); ++i)
            {
                const MipRecord &src = mipRecords[i];
                ReadMipRecord &dst = readMips[i];
                dst.mipLevel = src.mipLevel;
                dst.width = src.width;
                dst.height = src.height;
                dst.tileCountX = src.tileCountX;
                dst.tileCountY = src.tileCountY;
                dst.firstTileRecord = src.firstTileRecord;
                dst.tileRecordCount = src.tileRecordCount;
                dst.isTail = src.isTail;
                dst.byteOffset = src.byteOffset;
                dst.storedByteSize = src.byteSize;
                dst.rawByteSize = src.byteSize;
                dst.compression = kCompressionNone;
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
            for (size_t i = 0; i < tileRecords.size(); ++i)
            {
                const TileRecord &src = tileRecords[i];
                ReadTileRecord &dst = readTiles[i];
                dst.mipLevel = src.mipLevel;
                dst.tileX = src.tileX;
                dst.tileY = src.tileY;
                dst.width = src.width;
                dst.height = src.height;
                dst.byteOffset = src.byteOffset;
                dst.storedByteSize = src.byteSize;
                dst.rawByteSize = src.byteSize;
                dst.compression = kCompressionNone;
            }
        }
        else
        {
            std::vector<V4MipRecord> mipRecords(entry.mipRecordCount);
            in.seekg(static_cast<std::streamoff>(entry.mipRecordOffset), std::ios::beg);
            if (!detail::ReadBytes(&in, mipRecords.data(), mipRecords.size() * sizeof(V4MipRecord)))
            {
                if (outError)
                {
                    *outError = "failed reading v4 mip records: " + path.string();
                }
                return false;
            }
            for (size_t i = 0; i < mipRecords.size(); ++i)
            {
                const V4MipRecord &src = mipRecords[i];
                ReadMipRecord &dst = readMips[i];
                dst.mipLevel = src.mipLevel;
                dst.width = src.width;
                dst.height = src.height;
                dst.tileCountX = src.tileCountX;
                dst.tileCountY = src.tileCountY;
                dst.firstTileRecord = src.firstTileRecord;
                dst.tileRecordCount = src.tileRecordCount;
                dst.isTail = src.isTail;
                dst.byteOffset = src.byteOffset;
                dst.storedByteSize = src.storedByteSize;
                dst.rawByteSize = src.rawByteSize;
                dst.compression = src.compression;
            }
            std::vector<V4TileRecord> tileRecords(entry.tileRecordCount);
            in.seekg(static_cast<std::streamoff>(entry.tileRecordOffset), std::ios::beg);
            if (!detail::ReadBytes(&in, tileRecords.data(), tileRecords.size() * sizeof(V4TileRecord)))
            {
                if (outError)
                {
                    *outError = "failed reading v4 tile records: " + path.string();
                }
                return false;
            }
            for (size_t i = 0; i < tileRecords.size(); ++i)
            {
                const V4TileRecord &src = tileRecords[i];
                ReadTileRecord &dst = readTiles[i];
                dst.mipLevel = src.mipLevel;
                dst.tileX = src.tileX;
                dst.tileY = src.tileY;
                dst.width = src.width;
                dst.height = src.height;
                dst.byteOffset = src.byteOffset;
                dst.storedByteSize = src.storedByteSize;
                dst.rawByteSize = src.rawByteSize;
                dst.compression = src.compression;
            }
        }

        UdimImage image = {};
        if (!ReconstructBaseImageFromRecords(path, entry, readMips, readTiles, &in, &image, outError))
        {
            return false;
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
