#include "texture/virtual_texture/tile_file.h"

#include <cstring>

namespace ybi
{
namespace texture
{

namespace
{

static constexpr uint32_t kUdimMin = 1001;
static constexpr uint32_t kUdimMax = 1100;
static constexpr uint32_t kCompressionNone = 0u;

struct TileFileHeader
{
    char magic[8];
    uint32_t version = 3;
    uint32_t channels = 4;
    uint32_t elementType = 1;
    uint32_t udimCount = 0;
    uint64_t udimTableOffset = 0;
};

struct V3UdimEntry
{
    uint32_t udim = 0;
    uint32_t imageWidth = 0;
    uint32_t imageHeight = 0;
    uint32_t tileSize = 0;
    uint32_t mipCount = 0;
    uint32_t streamMipCount = 0;
    uint32_t tailMipCount = 0;
    uint32_t pixelFormat = 0;
    uint64_t mipRecordOffset = 0;
    uint32_t mipRecordCount = 0;
    uint32_t reserved1 = 0;
    uint64_t tileRecordOffset = 0;
    uint32_t tileRecordCount = 0;
    uint32_t reserved2 = 0;
    uint64_t payloadOffset = 0;
    uint64_t payloadBytes = 0;
};

struct V3MipRecord
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
    uint64_t byteSize = 0;
};

struct V3TileRecord
{
    uint32_t mipLevel = 0;
    uint32_t tileX = 0;
    uint32_t tileY = 0;
    uint32_t width = 0;
    uint32_t height = 0;
    uint64_t byteOffset = 0;
    uint64_t byteSize = 0;
};

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

bool ReadBytes(std::ifstream *stream, void *dst, size_t bytes)
{
    stream->read(reinterpret_cast<char *>(dst), static_cast<std::streamsize>(bytes));
    return stream->good();
}

bool ValidateBaseEntry(const V3UdimEntry &entry, const std::string &path, std::string *outError)
{
    if (entry.udim < kUdimMin || entry.udim > kUdimMax || entry.tileSize == 0 ||
        entry.mipCount == 0 || entry.mipRecordCount != entry.mipCount)
    {
        if (outError)
        {
            *outError = "invalid tile UDIM entry in: " + path;
        }
        return false;
    }
    return true;
}

TextureFormat ResolvePixelFormatForHeader(const TileFileHeader &header, const V3UdimEntry &entry)
{
    if (header.version >= 5 && IsValidTextureFormat(entry.pixelFormat))
    {
        return static_cast<TextureFormat>(entry.pixelFormat);
    }
    return TextureFormat::RGBA32_FLOAT;
}

bool OpenTileFileV3(const TileFileHeader &header,
                    VirtualTextureTileFile *outFile,
                    std::string *outError)
{
    std::vector<V3UdimEntry> entries(header.udimCount);
    outFile->stream.seekg(static_cast<std::streamoff>(header.udimTableOffset), std::ios::beg);
    if (!ReadBytes(&outFile->stream, entries.data(), entries.size() * sizeof(V3UdimEntry)))
    {
        if (outError)
        {
            *outError = "failed reading tile file UDIM table: " + outFile->path;
        }
        return false;
    }

    for (const V3UdimEntry &entry : entries)
    {
        if (!ValidateBaseEntry(entry, outFile->path, outError))
        {
            return false;
        }

        std::vector<V3MipRecord> mipRecords(entry.mipRecordCount);
        outFile->stream.seekg(static_cast<std::streamoff>(entry.mipRecordOffset), std::ios::beg);
        if (!ReadBytes(&outFile->stream, mipRecords.data(), mipRecords.size() * sizeof(V3MipRecord)))
        {
            if (outError)
            {
                *outError = "failed reading mip records in: " + outFile->path;
            }
            return false;
        }

        std::vector<V3TileRecord> tileRecords(entry.tileRecordCount);
        outFile->stream.seekg(static_cast<std::streamoff>(entry.tileRecordOffset), std::ios::beg);
        if (!ReadBytes(&outFile->stream, tileRecords.data(), tileRecords.size() * sizeof(V3TileRecord)))
        {
            if (outError)
            {
                *outError = "failed reading tile records in: " + outFile->path;
            }
            return false;
        }

        VirtualTextureUdimTable table = {};
        table.imageWidth = entry.imageWidth;
        table.imageHeight = entry.imageHeight;
        table.tileSize = entry.tileSize;
        table.pixelFormat = ResolvePixelFormatForHeader(header, entry);
        table.mips.resize(entry.mipCount);

        for (const V3MipRecord &srcMip : mipRecords)
        {
            if (srcMip.mipLevel >= entry.mipCount)
            {
                if (outError)
                {
                    *outError = "mip level out of range in: " + outFile->path;
                }
                return false;
            }

            VirtualTextureMipTable &dstMip = table.mips[srcMip.mipLevel];
            if (dstMip.width != 0 || dstMip.height != 0)
            {
                if (outError)
                {
                    *outError = "duplicate mip level in: " + outFile->path;
                }
                return false;
            }

            dstMip.level = srcMip.mipLevel;
            dstMip.width = srcMip.width;
            dstMip.height = srcMip.height;
            dstMip.tileCountX = srcMip.tileCountX;
            dstMip.tileCountY = srcMip.tileCountY;
            dstMip.isTail = (srcMip.isTail != 0u);
            dstMip.tailByteOffset = srcMip.byteOffset;
            dstMip.tailStoredByteSize = srcMip.byteSize;
            dstMip.tailRawByteSize = srcMip.byteSize;
            dstMip.tailCompression = kCompressionNone;

            if (dstMip.isTail)
            {
                const uint64_t expectedBytes = static_cast<uint64_t>(dstMip.width) *
                                               static_cast<uint64_t>(dstMip.height) * 4u * sizeof(float);
                if (dstMip.tailRawByteSize != expectedBytes)
                {
                    if (outError)
                    {
                        *outError = "tail mip byte size mismatch in: " + outFile->path;
                    }
                    return false;
                }
                continue;
            }

            if (srcMip.firstTileRecord + srcMip.tileRecordCount > tileRecords.size())
            {
                if (outError)
                {
                    *outError = "tile record range out of bounds in: " + outFile->path;
                }
                return false;
            }

            dstMip.records.resize(srcMip.tileRecordCount);
            for (uint32_t i = 0u; i < srcMip.tileRecordCount; ++i)
            {
                const V3TileRecord &srcTile = tileRecords[srcMip.firstTileRecord + i];
                VirtualTextureTileRecord &dstTile = dstMip.records[i];
                dstTile.mipLevel = srcTile.mipLevel;
                dstTile.tileX = srcTile.tileX;
                dstTile.tileY = srcTile.tileY;
                dstTile.width = srcTile.width;
                dstTile.height = srcTile.height;
                dstTile.byteOffset = srcTile.byteOffset;
                dstTile.storedByteSize = srcTile.byteSize;
                dstTile.rawByteSize = srcTile.byteSize;
                dstTile.compression = kCompressionNone;
            }
        }

        for (size_t mip = 0; mip < table.mips.size(); ++mip)
        {
            if (table.mips[mip].width == 0 || table.mips[mip].height == 0)
            {
                if (outError)
                {
                    *outError = "missing mip level in: " + outFile->path;
                }
                return false;
            }
        }

        outFile->totalTextureBytes += static_cast<uint64_t>(entry.imageWidth) *
                                      static_cast<uint64_t>(entry.imageHeight) *
                                      static_cast<uint64_t>(
                                          TextureFormatPixelBytes(table.pixelFormat));
        outFile->udims.emplace(entry.udim, std::move(table));
    }

    return true;
}

bool OpenTileFileV4(const TileFileHeader &header,
                    VirtualTextureTileFile *outFile,
                    std::string *outError)
{
    std::vector<V3UdimEntry> entries(header.udimCount);
    outFile->stream.seekg(static_cast<std::streamoff>(header.udimTableOffset), std::ios::beg);
    if (!ReadBytes(&outFile->stream, entries.data(), entries.size() * sizeof(V3UdimEntry)))
    {
        if (outError)
        {
            *outError = "failed reading tile file v4 UDIM table: " + outFile->path;
        }
        return false;
    }

    for (const V3UdimEntry &entry : entries)
    {
        if (!ValidateBaseEntry(entry, outFile->path, outError))
        {
            return false;
        }

        std::vector<V4MipRecord> mipRecords(entry.mipRecordCount);
        outFile->stream.seekg(static_cast<std::streamoff>(entry.mipRecordOffset), std::ios::beg);
        if (!ReadBytes(&outFile->stream, mipRecords.data(), mipRecords.size() * sizeof(V4MipRecord)))
        {
            if (outError)
            {
                *outError = "failed reading v4 mip records in: " + outFile->path;
            }
            return false;
        }

        std::vector<V4TileRecord> tileRecords(entry.tileRecordCount);
        outFile->stream.seekg(static_cast<std::streamoff>(entry.tileRecordOffset), std::ios::beg);
        if (!ReadBytes(&outFile->stream, tileRecords.data(), tileRecords.size() * sizeof(V4TileRecord)))
        {
            if (outError)
            {
                *outError = "failed reading v4 tile records in: " + outFile->path;
            }
            return false;
        }

        VirtualTextureUdimTable table = {};
        table.imageWidth = entry.imageWidth;
        table.imageHeight = entry.imageHeight;
        table.tileSize = entry.tileSize;
        table.pixelFormat = ResolvePixelFormatForHeader(header, entry);
        table.mips.resize(entry.mipCount);

        for (const V4MipRecord &srcMip : mipRecords)
        {
            if (srcMip.mipLevel >= entry.mipCount)
            {
                if (outError)
                {
                    *outError = "mip level out of range in: " + outFile->path;
                }
                return false;
            }

            VirtualTextureMipTable &dstMip = table.mips[srcMip.mipLevel];
            if (dstMip.width != 0 || dstMip.height != 0)
            {
                if (outError)
                {
                    *outError = "duplicate mip level in: " + outFile->path;
                }
                return false;
            }

            dstMip.level = srcMip.mipLevel;
            dstMip.width = srcMip.width;
            dstMip.height = srcMip.height;
            dstMip.tileCountX = srcMip.tileCountX;
            dstMip.tileCountY = srcMip.tileCountY;
            dstMip.isTail = (srcMip.isTail != 0u);
            dstMip.tailByteOffset = srcMip.byteOffset;
            dstMip.tailStoredByteSize = srcMip.storedByteSize;
            dstMip.tailRawByteSize = srcMip.rawByteSize;
            dstMip.tailCompression = srcMip.compression;

            if (dstMip.isTail)
            {
                const uint64_t expectedBytes = static_cast<uint64_t>(dstMip.width) *
                                               static_cast<uint64_t>(dstMip.height) * 4u * sizeof(float);
                if (dstMip.tailRawByteSize != expectedBytes)
                {
                    if (outError)
                    {
                        *outError = "tail mip raw byte size mismatch in: " + outFile->path;
                    }
                    return false;
                }
                continue;
            }

            if (srcMip.firstTileRecord + srcMip.tileRecordCount > tileRecords.size())
            {
                if (outError)
                {
                    *outError = "tile record range out of bounds in: " + outFile->path;
                }
                return false;
            }

            dstMip.records.resize(srcMip.tileRecordCount);
            for (uint32_t i = 0u; i < srcMip.tileRecordCount; ++i)
            {
                const V4TileRecord &srcTile = tileRecords[srcMip.firstTileRecord + i];
                VirtualTextureTileRecord &dstTile = dstMip.records[i];
                dstTile.mipLevel = srcTile.mipLevel;
                dstTile.tileX = srcTile.tileX;
                dstTile.tileY = srcTile.tileY;
                dstTile.width = srcTile.width;
                dstTile.height = srcTile.height;
                dstTile.byteOffset = srcTile.byteOffset;
                dstTile.storedByteSize = srcTile.storedByteSize;
                dstTile.rawByteSize = srcTile.rawByteSize;
                dstTile.compression = srcTile.compression;
            }
        }

        for (size_t mip = 0; mip < table.mips.size(); ++mip)
        {
            if (table.mips[mip].width == 0 || table.mips[mip].height == 0)
            {
                if (outError)
                {
                    *outError = "missing mip level in: " + outFile->path;
                }
                return false;
            }
        }

        outFile->totalTextureBytes += static_cast<uint64_t>(entry.imageWidth) *
                                      static_cast<uint64_t>(entry.imageHeight) *
                                      static_cast<uint64_t>(
                                          TextureFormatPixelBytes(table.pixelFormat));
        outFile->udims.emplace(entry.udim, std::move(table));
    }

    return true;
}

} // namespace

bool OpenVirtualTextureTileFile(const std::string &path,
                                VirtualTextureTileFile *outFile,
                                std::string *outError)
{
    if (!outFile)
    {
        if (outError)
        {
            *outError = "null output tile file handle";
        }
        return false;
    }

    outFile->path = path;
    outFile->stream = std::ifstream(path, std::ios::binary);
    outFile->udims.clear();
    outFile->totalTextureBytes = 0;
    if (!outFile->stream.is_open())
    {
        if (outError)
        {
            *outError = "failed opening tile file: " + path;
        }
        return false;
    }

    TileFileHeader header = {};
    if (!ReadBytes(&outFile->stream, &header, sizeof(header)))
    {
        if (outError)
        {
            *outError = "failed reading tile file header: " + path;
        }
        return false;
    }

    if (header.channels != 4 || header.elementType != 1 || header.udimCount == 0)
    {
        if (outError)
        {
            *outError = "invalid tile file header: " + path;
        }
        return false;
    }

    bool opened = false;
    if (std::memcmp(header.magic, "YBITILE5", 8) == 0 && header.version == 5)
    {
        opened = OpenTileFileV4(header, outFile, outError);
    }
    else if (std::memcmp(header.magic, "YBITILE4", 8) == 0 && header.version == 4)
    {
        opened = OpenTileFileV4(header, outFile, outError);
    }
    else if (std::memcmp(header.magic, "YBITILE3", 8) == 0 && header.version == 3)
    {
        opened = OpenTileFileV3(header, outFile, outError);
    }
    else if (outError)
    {
        *outError = "unsupported tile file version: " + path;
    }

    outFile->stream.close();
    return opened;
}

} // namespace texture
} // namespace ybi
