#pragma once

#include "exr_mips.h"
#include "tile_binary.h"

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace ybi
{
namespace tileprep
{

enum class TextureSemanticClass : uint8_t
{
    Unknown = 0,
    Scalar = 1,
    Normal = 2,
    Color = 3,
};

bool ChoosePixelFormat(TextureSemanticClass semanticClass,
                       uint32_t sourceChannelCount,
                       ExrNumericType sourceNumericType,
                       texture::TextureFormat *outFormat,
                       std::string *outError);

bool ConvertMipChainToPixelFormat(texture::TextureFormat pixelFormat,
                                  std::vector<tilebin::UdimMipImage> *mipLevels,
                                  std::string *outError);

bool ConvertPixelsToPixelFormat(texture::TextureFormat pixelFormat,
                                const std::vector<float> &rgbaPixels,
                                std::vector<float> *outPixels,
                                std::string *outError);

void ExpandPixelsToRgba(texture::TextureFormat pixelFormat,
                        const std::vector<float> &pixels,
                        std::vector<float> *outRgba);

bool EncodePixelPayload(texture::TextureFormat pixelFormat,
                        const std::vector<float> &pixels,
                        std::vector<unsigned char> *outBytes,
                        std::string *outError);

bool DecodePixelPayload(texture::TextureFormat pixelFormat,
                        const unsigned char *bytes,
                        size_t byteCount,
                        std::vector<float> *outPixels,
                        std::string *outError);

} // namespace tileprep
} // namespace ybi
