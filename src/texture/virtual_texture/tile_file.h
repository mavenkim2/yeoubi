#pragma once

#include <cstdint>
#include <fstream>
#include <string>
#include <unordered_map>
#include <vector>

namespace ybi
{
namespace texture
{

struct VirtualTextureTileRecord
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
    uint32_t reserved0 = 0;
};

struct VirtualTextureMipTable
{
    uint32_t level = 0;
    uint32_t width = 0;
    uint32_t height = 0;
    uint32_t tileCountX = 0;
    uint32_t tileCountY = 0;
    bool isTail = false;
    uint64_t tailByteOffset = 0;
    uint64_t tailStoredByteSize = 0;
    uint64_t tailRawByteSize = 0;
    uint32_t tailCompression = 0;
    std::vector<VirtualTextureTileRecord> records;
};

struct VirtualTextureUdimTable
{
    uint32_t imageWidth = 0;
    uint32_t imageHeight = 0;
    uint32_t tileSize = 0;
    std::vector<VirtualTextureMipTable> mips;
};

struct VirtualTextureTileFile
{
    std::string path;
    std::ifstream stream;
    uint64_t totalTextureBytes = 0;
    std::unordered_map<uint32_t, VirtualTextureUdimTable> udims;
};

bool OpenVirtualTextureTileFile(const std::string &path,
                                VirtualTextureTileFile *outFile,
                                std::string *outError);

bool ReadVirtualTextureTile(VirtualTextureTileFile *file,
                            uint32_t udim,
                            uint32_t mip,
                            uint32_t tileX,
                            uint32_t tileY,
                            std::vector<unsigned char> *outRgba8,
                            uint32_t *outWidth,
                            uint32_t *outHeight,
                            uint64_t *outSourceBytes,
                            std::string *outError);

bool ReadVirtualTextureTailMip(VirtualTextureTileFile *file,
                               uint32_t udim,
                               uint32_t maxDim,
                               std::vector<unsigned char> *outRgba8,
                               uint32_t *outWidth,
                               uint32_t *outHeight,
                               uint32_t *outMip,
                               uint64_t *outSourceBytes,
                               std::string *outError);

} // namespace texture
} // namespace ybi
