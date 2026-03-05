#include "texture/virtual_texture/tile_file.h"

#include <algorithm>
#include <cmath>

namespace ybi
{
namespace texture
{

namespace
{

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
        const uint64_t expectedBytes = static_cast<uint64_t>(mipTable.width) *
                                       static_cast<uint64_t>(mipTable.height) * 4u * sizeof(float);
        if (mipTable.tailByteSize != expectedBytes)
        {
            if (outError)
            {
                *outError = "tail mip byte size mismatch";
            }
            return false;
        }

        std::vector<float> rgbaF32(static_cast<size_t>(mipTable.width) *
                                   static_cast<size_t>(mipTable.height) * 4u);
        file->stream.clear();
        file->stream.seekg(static_cast<std::streamoff>(mipTable.tailByteOffset), std::ios::beg);
        if (!ReadBytes(&file->stream, rgbaF32.data(), static_cast<size_t>(mipTable.tailByteSize)))
        {
            if (outError)
            {
                *outError = "failed reading tail mip payload";
            }
            return false;
        }

        ConvertRgbaF32ToU8(rgbaF32, outRgba8);
        *outWidth = mipTable.width;
        *outHeight = mipTable.height;
        *outSourceBytes = mipTable.tailByteSize;
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

    ConvertRgbaF32ToU8(rgbaF32, outRgba8);
    *outWidth = record.width;
    *outHeight = record.height;
    *outSourceBytes = record.byteSize;
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
