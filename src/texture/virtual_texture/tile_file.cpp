#include "texture/virtual_texture/tile_file.h"

#include "miniz.h"
#include "util/half_float.h"

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

bool ReadDecodedPayload(VirtualTextureTileFile *file,
                        uint64_t byteOffset,
                        uint64_t storedByteSize,
                        uint64_t rawByteSize,
                        uint32_t compression,
                        std::vector<unsigned char> *outRaw,
                        std::string *outError)
{
    if (!EnsureMzUlongRange(storedByteSize, outError) || !EnsureMzUlongRange(rawByteSize, outError))
    {
        return false;
    }

    std::ifstream stream(file->path, std::ios::binary);
    if (!stream.is_open())
    {
        if (outError)
        {
            *outError = "failed opening tile file for payload read: " + file->path;
        }
        return false;
    }

    std::vector<unsigned char> stored(static_cast<size_t>(storedByteSize));
    stream.seekg(static_cast<std::streamoff>(byteOffset), std::ios::beg);
    if (!ReadBytes(&stream, stored.data(), stored.size()))
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
        *outRaw = std::move(stored);
        return true;
    }
    else if (compression == kCompressionDeflate)
    {
        outRaw->assign(static_cast<size_t>(rawByteSize), 0u);
        mz_ulong decodedBytes = static_cast<mz_ulong>(rawByteSize);
        const int rc = mz_uncompress(outRaw->data(),
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
    return true;
}

bool DecodePixelsToFloats(VirtualTexturePixelFormat pixelFormat,
                          const std::vector<unsigned char> &rawBytes,
                          std::vector<float> *outPixels,
                          std::string *outError)
{
    const uint32_t channelCount = VirtualTexturePixelFormatChannelCount(pixelFormat);
    const uint32_t bytesPerChannel = VirtualTexturePixelFormatBytesPerChannel(pixelFormat);
    if (channelCount == 0u || bytesPerChannel == 0u)
    {
        if (outError)
        {
            *outError = "invalid virtual texture pixel format";
        }
        return false;
    }
    if ((rawBytes.size() % static_cast<size_t>(bytesPerChannel)) != 0u)
    {
        if (outError)
        {
            *outError = "typed tile payload size mismatch";
        }
        return false;
    }

    const size_t sampleCount = rawBytes.size() / static_cast<size_t>(bytesPerChannel);
    outPixels->assign(sampleCount, 0.0f);
    if (bytesPerChannel == sizeof(float))
    {
        std::memcpy(outPixels->data(), rawBytes.data(), rawBytes.size());
        return true;
    }

    for (size_t i = 0; i < sampleCount; ++i)
    {
        uint16_t halfBits = 0u;
        std::memcpy(&halfBits,
                    rawBytes.data() + i * sizeof(uint16_t),
                    sizeof(uint16_t));
        (*outPixels)[i] = util::HalfBitsToFloat(halfBits);
    }
    return true;
}

void ExpandTypedPixelsToRgba8(VirtualTexturePixelFormat pixelFormat,
                              const std::vector<float> &pixels,
                              std::vector<unsigned char> *outRgba8)
{
    const uint32_t channelCount = VirtualTexturePixelFormatChannelCount(pixelFormat);
    if (channelCount == 0u || (pixels.size() % static_cast<size_t>(channelCount)) != 0u)
    {
        outRgba8->clear();
        return;
    }

    const size_t pixelCount = pixels.size() / static_cast<size_t>(channelCount);
    outRgba8->assign(pixelCount * 4u, 0u);
    for (size_t pixelIndex = 0; pixelIndex < pixelCount; ++pixelIndex)
    {
        const float *src = pixels.data() + pixelIndex * static_cast<size_t>(channelCount);
        float rgba[4] = {0.0f, 0.0f, 0.0f, 1.0f};
        switch (channelCount)
        {
            case 1u:
                rgba[0] = src[0];
                rgba[1] = src[0];
                rgba[2] = src[0];
                break;
            case 2u:
                rgba[0] = src[0];
                rgba[1] = src[1];
                break;
            case 4u:
                std::memcpy(rgba, src, sizeof(rgba));
                break;
            default:
                outRgba8->clear();
                return;
        }

        for (uint32_t c = 0u; c < 4u; ++c)
        {
            const float v = std::min(1.0f, std::max(0.0f, rgba[c]));
            (*outRgba8)[pixelIndex * 4u + c] =
                static_cast<unsigned char>(std::round(v * 255.0f));
        }
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
    const uint32_t channelCount = VirtualTexturePixelFormatChannelCount(table.pixelFormat);
    const uint32_t bytesPerChannel = VirtualTexturePixelFormatBytesPerChannel(table.pixelFormat);
    if (channelCount == 0u || bytesPerChannel == 0u)
    {
        if (outError)
        {
            *outError = "invalid typed virtual texture format";
        }
        return false;
    }

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
                                          static_cast<uint64_t>(mipTable.height) *
                                          static_cast<uint64_t>(channelCount) *
                                          static_cast<uint64_t>(bytesPerChannel);
        if (mipTable.tailRawByteSize != expectedRawBytes)
        {
            if (outError)
            {
                *outError = "tail mip byte size mismatch";
            }
            return false;
        }

        std::vector<unsigned char> rawPayload;
        if (!ReadDecodedPayload(file,
                                mipTable.tailByteOffset,
                                mipTable.tailStoredByteSize,
                                mipTable.tailRawByteSize,
                                mipTable.tailCompression,
                                &rawPayload,
                                outError))
        {
            return false;
        }

        std::vector<float> pixels;
        if (!DecodePixelsToFloats(table.pixelFormat, rawPayload, &pixels, outError))
        {
            return false;
        }
        ExpandTypedPixelsToRgba8(table.pixelFormat, pixels, outRgba8);
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
                                      static_cast<uint64_t>(record.height) *
                                      static_cast<uint64_t>(channelCount) *
                                      static_cast<uint64_t>(bytesPerChannel);
    if (record.rawByteSize != expectedRawBytes)
    {
        if (outError)
        {
            *outError = "tile byte size mismatch";
        }
        return false;
    }

    std::vector<unsigned char> rawPayload;
    if (!ReadDecodedPayload(file,
                            record.byteOffset,
                            record.storedByteSize,
                            record.rawByteSize,
                            record.compression,
                            &rawPayload,
                            outError))
    {
        return false;
    }

    std::vector<float> pixels;
    if (!DecodePixelsToFloats(table.pixelFormat, rawPayload, &pixels, outError))
    {
        return false;
    }
    ExpandTypedPixelsToRgba8(table.pixelFormat, pixels, outRgba8);
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
