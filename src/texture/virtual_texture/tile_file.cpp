#include "texture/virtual_texture/tile_file.h"

#include "miniz.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>

namespace ybi
{
namespace texture
{

namespace
{

static constexpr uint32_t kCompressionNone = 0u;
static constexpr uint32_t kCompressionDeflate = 1u;

bool ReadBytes(std::ifstream *stream, void *dst, size_t bytes)
{
    stream->read(reinterpret_cast<char *>(dst), static_cast<std::streamsize>(bytes));
    return stream->good();
}

void ConvertRgbaF32ToU8(const std::vector<float> &rgbaF32, std::vector<unsigned char> *outRgba8)
{
    outRgba8->resize(rgbaF32.size());
    for (size_t i = 0; i < rgbaF32.size(); ++i)
    {
        const float v = std::min(1.0f, std::max(0.0f, rgbaF32[i]));
        (*outRgba8)[i] = static_cast<unsigned char>(std::round(v * 255.0f));
    }
}

bool EnsureMzUlongRange(uint64_t bytes, std::string *outError)
{
    if (bytes > static_cast<uint64_t>(std::numeric_limits<mz_ulong>::max()))
    {
        if (outError)
        {
            *outError = "payload too large for miniz";
        }
        return false;
    }
    return true;
}

bool ReadPayloadAsFloats(VirtualTextureTileFile *file,
                         uint64_t byteOffset,
                         uint64_t storedByteSize,
                         uint64_t rawByteSize,
                         uint32_t compression,
                         std::vector<float> *outRgbaF32,
                         std::string *outError)
{
    if ((rawByteSize % sizeof(float)) != 0u)
    {
        if (outError)
        {
            *outError = "raw payload byte size is not float-aligned";
        }
        return false;
    }
    if (!EnsureMzUlongRange(storedByteSize, outError) || !EnsureMzUlongRange(rawByteSize, outError))
    {
        return false;
    }

    std::vector<unsigned char> stored(static_cast<size_t>(storedByteSize));
    file->stream.clear();
    file->stream.seekg(static_cast<std::streamoff>(byteOffset), std::ios::beg);
    if (!ReadBytes(&file->stream, stored.data(), stored.size()))
    {
        if (outError)
        {
            *outError = "failed reading tile payload";
        }
        return false;
    }

    std::vector<unsigned char> raw;
    if (compression == kCompressionNone)
    {
        if (storedByteSize != rawByteSize)
        {
            if (outError)
            {
                *outError = "raw payload size mismatch";
            }
            return false;
        }
        raw = std::move(stored);
    }
    else if (compression == kCompressionDeflate)
    {
        raw.assign(static_cast<size_t>(rawByteSize), 0u);
        mz_ulong decodedBytes = static_cast<mz_ulong>(rawByteSize);
        const int rc = mz_uncompress(raw.data(),
                                     &decodedBytes,
                                     stored.data(),
                                     static_cast<mz_ulong>(storedByteSize));
        if (rc != MZ_OK || decodedBytes != static_cast<mz_ulong>(rawByteSize))
        {
            if (outError)
            {
                *outError = "failed decompressing tile payload";
            }
            return false;
        }
    }
    else
    {
        if (outError)
        {
            *outError = "unsupported tile payload compression";
        }
        return false;
    }

    outRgbaF32->assign(rawByteSize / sizeof(float), 0.0f);
    std::memcpy(outRgbaF32->data(), raw.data(), static_cast<size_t>(rawByteSize));
    return true;
}

} // namespace

bool ReadVirtualTextureTile(VirtualTextureTileFile *file,
                            uint32_t udim,
                            uint32_t mip,
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
    if (mip >= table.mips.size())
    {
        if (outError)
        {
            *outError = "mip level out of range";
        }
        return false;
    }
    const VirtualTextureMipTable &mipTable = table.mips[mip];

    if (mipTable.isTail)
    {
        if (tileX != 0u || tileY != 0u)
        {
            if (outError)
            {
                *outError = "tail mip only supports tile (0,0)";
            }
            return false;
        }
        const uint64_t expectedRawBytes = static_cast<uint64_t>(mipTable.width) *
                                          static_cast<uint64_t>(mipTable.height) * 4u * sizeof(float);
        if (mipTable.tailRawByteSize != expectedRawBytes)
        {
            if (outError)
            {
                *outError = "tail mip byte size mismatch";
            }
            return false;
        }

        std::vector<float> rgbaF32;
        if (!ReadPayloadAsFloats(file,
                                 mipTable.tailByteOffset,
                                 mipTable.tailStoredByteSize,
                                 mipTable.tailRawByteSize,
                                 mipTable.tailCompression,
                                 &rgbaF32,
                                 outError))
        {
            return false;
        }

        ConvertRgbaF32ToU8(rgbaF32, outRgba8);
        *outWidth = mipTable.width;
        *outHeight = mipTable.height;
        *outSourceBytes = mipTable.tailStoredByteSize;
        return true;
    }

    if (tileX >= mipTable.tileCountX || tileY >= mipTable.tileCountY)
    {
        if (outError)
        {
            *outError = "tile coordinates out of range";
        }
        return false;
    }
    const uint64_t tileIndex = static_cast<uint64_t>(tileY) * static_cast<uint64_t>(mipTable.tileCountX) +
                               static_cast<uint64_t>(tileX);
    if (tileIndex >= mipTable.records.size())
    {
        if (outError)
        {
            *outError = "tile index out of range";
        }
        return false;
    }

    const VirtualTextureTileRecord &record = mipTable.records[static_cast<size_t>(tileIndex)];
    if (record.width == 0 || record.height == 0)
    {
        if (outError)
        {
            *outError = "tile record has zero dimension";
        }
        return false;
    }

    const uint64_t expectedRawBytes = static_cast<uint64_t>(record.width) *
                                      static_cast<uint64_t>(record.height) * 4u * sizeof(float);
    if (record.rawByteSize != expectedRawBytes)
    {
        if (outError)
        {
            *outError = "tile byte size mismatch";
        }
        return false;
    }

    std::vector<float> rgbaF32;
    if (!ReadPayloadAsFloats(file,
                             record.byteOffset,
                             record.storedByteSize,
                             record.rawByteSize,
                             record.compression,
                             &rgbaF32,
                             outError))
    {
        return false;
    }

    ConvertRgbaF32ToU8(rgbaF32, outRgba8);
    *outWidth = record.width;
    *outHeight = record.height;
    *outSourceBytes = record.storedByteSize;
    return true;
}

bool ReadVirtualTextureTailMip(VirtualTextureTileFile *file,
                               uint32_t udim,
                               uint32_t maxDim,
                               std::vector<unsigned char> *outRgba8,
                               uint32_t *outWidth,
                               uint32_t *outHeight,
                               uint32_t *outMip,
                               uint64_t *outSourceBytes,
                               std::string *outError)
{
    if (!file || !outRgba8 || !outWidth || !outHeight || !outMip || !outSourceBytes)
    {
        if (outError)
        {
            *outError = "invalid read tail arguments";
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
    const VirtualTextureMipTable *bestWithin = nullptr;
    const VirtualTextureMipTable *bestFallback = nullptr;
    for (const VirtualTextureMipTable &mipTable : table.mips)
    {
        if (!mipTable.isTail)
        {
            continue;
        }
        const uint32_t dim = std::max(mipTable.width, mipTable.height);
        if (!bestFallback || dim > std::max(bestFallback->width, bestFallback->height))
        {
            bestFallback = &mipTable;
        }
        if (dim <= maxDim)
        {
            if (!bestWithin || dim > std::max(bestWithin->width, bestWithin->height))
            {
                bestWithin = &mipTable;
            }
        }
    }

    const VirtualTextureMipTable *chosen = bestWithin ? bestWithin : bestFallback;
    if (!chosen)
    {
        if (outError)
        {
            *outError = "no tail mip present for udim " + std::to_string(udim);
        }
        return false;
    }

    if (!ReadVirtualTextureTile(file,
                                udim,
                                chosen->level,
                                0u,
                                0u,
                                outRgba8,
                                outWidth,
                                outHeight,
                                outSourceBytes,
                                outError))
    {
        return false;
    }

    *outMip = chosen->level;
    return true;
}

} // namespace texture
} // namespace ybi
