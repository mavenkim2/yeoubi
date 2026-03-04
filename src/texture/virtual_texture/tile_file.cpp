#include "texture/virtual_texture/tile_file.h"

#include <algorithm>
#include <cmath>
#include <cstring>

namespace ybi
{
namespace texture
{

namespace
{
static constexpr uint32_t kUdimMin = 1001;
static constexpr uint32_t kUdimMax = 1100;

struct TileFileHeader
{
    char magic[8];
    uint32_t version = 2;
    uint32_t channels = 4;
    uint32_t elementType = 1;
    uint32_t udimCount = 0;
    uint64_t udimTableOffset = 0;
};

struct UdimEntry
{
    uint32_t udim = 0;
    uint32_t imageWidth = 0;
    uint32_t imageHeight = 0;
    uint32_t tileSize = 0;
    uint32_t tileCountX = 0;
    uint32_t tileCountY = 0;
    uint32_t tileCount = 0;
    uint64_t tileRecordOffset = 0;
    uint32_t tileRecordCount = 0;
    uint32_t reserved0 = 0;
    uint64_t payloadOffset = 0;
    uint64_t payloadBytes = 0;
};

struct TileRecord
{
    uint32_t tileX = 0;
    uint32_t tileY = 0;
    uint32_t width = 0;
    uint32_t height = 0;
    uint64_t byteOffset = 0;
    uint64_t byteSize = 0;
};

bool ReadBytes(std::ifstream *stream, void *dst, size_t bytes)
{
    stream->read(reinterpret_cast<char *>(dst), static_cast<std::streamsize>(bytes));
    return stream->good();
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

    if (std::memcmp(header.magic, "YBITILE2", 8) != 0 || header.version != 2 ||
        header.channels != 4 || header.elementType != 1 || header.udimCount == 0)
    {
        if (outError)
        {
            *outError = "invalid tile file header: " + path;
        }
        return false;
    }

    std::vector<UdimEntry> entries(header.udimCount);
    outFile->stream.seekg(static_cast<std::streamoff>(header.udimTableOffset), std::ios::beg);
    if (!ReadBytes(&outFile->stream, entries.data(), entries.size() * sizeof(UdimEntry)))
    {
        if (outError)
        {
            *outError = "failed reading tile file UDIM table: " + path;
        }
        return false;
    }

    for (const UdimEntry &entry : entries)
    {
        if (entry.udim < kUdimMin || entry.udim > kUdimMax || entry.tileSize == 0 ||
            entry.tileCount == 0 || entry.tileRecordCount != entry.tileCount ||
            entry.tileCount != entry.tileCountX * entry.tileCountY)
        {
            if (outError)
            {
                *outError = "invalid tile UDIM entry in: " + path;
            }
            return false;
        }

        std::vector<TileRecord> records(entry.tileRecordCount);
        outFile->stream.seekg(static_cast<std::streamoff>(entry.tileRecordOffset), std::ios::beg);
        if (!ReadBytes(&outFile->stream, records.data(), records.size() * sizeof(TileRecord)))
        {
            if (outError)
            {
                *outError = "failed reading tile records in: " + path;
            }
            return false;
        }

        VirtualTextureUdimTable table = {};
        table.imageWidth = entry.imageWidth;
        table.imageHeight = entry.imageHeight;
        table.tileSize = entry.tileSize;
        table.tileCountX = entry.tileCountX;
        table.tileCountY = entry.tileCountY;
        table.records.resize(records.size());
        for (size_t i = 0; i < records.size(); ++i)
        {
            const TileRecord &src = records[i];
            VirtualTextureTileRecord &dst = table.records[i];
            dst.width = src.width;
            dst.height = src.height;
            dst.byteOffset = src.byteOffset;
            dst.byteSize = src.byteSize;
        }
        outFile->totalTextureBytes += static_cast<uint64_t>(entry.imageWidth) *
                                      static_cast<uint64_t>(entry.imageHeight) * 4u * sizeof(float);
        outFile->udims.emplace(entry.udim, std::move(table));
    }

    return true;
}

bool ReadVirtualTextureTile(VirtualTextureTileFile *file,
                            uint32_t udim,
                            uint32_t tileX,
                            uint32_t tileY,
                            std::vector<unsigned char> *outRgba8,
                            uint32_t *outWidth,
                            uint32_t *outHeight,
                            uint64_t *outSourceBytes,
                            std::string *outError)
{
    if (!file || !outRgba8 || !outWidth || !outHeight || !outSourceBytes)
    {
        if (outError)
        {
            *outError = "invalid read tile arguments";
        }
        return false;
    }

    auto udimIter = file->udims.find(udim);
    if (udimIter == file->udims.end())
    {
        if (outError)
        {
            *outError = "udim not found in tile file: " + std::to_string(udim);
        }
        return false;
    }

    const VirtualTextureUdimTable &table = udimIter->second;
    if (tileX >= table.tileCountX || tileY >= table.tileCountY)
    {
        if (outError)
        {
            *outError = "tile coordinates out of range";
        }
        return false;
    }
    const uint64_t tileIndex = static_cast<uint64_t>(tileY) * static_cast<uint64_t>(table.tileCountX) +
                               static_cast<uint64_t>(tileX);
    if (tileIndex >= table.records.size())
    {
        if (outError)
        {
            *outError = "tile index out of range";
        }
        return false;
    }

    const VirtualTextureTileRecord &record = table.records[static_cast<size_t>(tileIndex)];
    if (record.width == 0 || record.height == 0)
    {
        if (outError)
        {
            *outError = "tile record has zero dimension";
        }
        return false;
    }

    const uint64_t expectedBytes = static_cast<uint64_t>(record.width) *
                                   static_cast<uint64_t>(record.height) * 4u * sizeof(float);
    if (record.byteSize != expectedBytes)
    {
        if (outError)
        {
            *outError = "tile byte size mismatch";
        }
        return false;
    }

    std::vector<float> rgbaF32(static_cast<size_t>(record.width) *
                               static_cast<size_t>(record.height) * 4u);
    file->stream.clear();
    file->stream.seekg(static_cast<std::streamoff>(record.byteOffset), std::ios::beg);
    if (!ReadBytes(&file->stream, rgbaF32.data(), static_cast<size_t>(record.byteSize)))
    {
        if (outError)
        {
            *outError = "failed reading tile payload";
        }
        return false;
    }

    outRgba8->resize(rgbaF32.size());
    for (size_t i = 0; i < rgbaF32.size(); ++i)
    {
        const float v = std::min(1.0f, std::max(0.0f, rgbaF32[i]));
        (*outRgba8)[i] = static_cast<unsigned char>(std::round(v * 255.0f));
    }
    *outWidth = record.width;
    *outHeight = record.height;
    *outSourceBytes = record.byteSize;
    return true;
}

} // namespace texture
} // namespace ybi
