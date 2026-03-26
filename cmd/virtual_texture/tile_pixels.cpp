#include "tile_pixels.h"

#include "util/half_float.h"

#include <cstring>

namespace ybi
{
namespace tileprep
{

namespace
{

uint32_t ChannelCountForFormat(TextureFormat pixelFormat)
{
    return TextureFormatChannelCount(pixelFormat);
}

bool ValidateRgbaPixels(const std::vector<float> &rgbaPixels, std::string *outError)
{
    if ((rgbaPixels.size() % 4u) != 0u)
    {
        if (outError)
        {
            *outError = "pixel buffer is not packed RGBA float data";
        }
        return false;
    }
    return true;
}

} // namespace

bool ChoosePixelFormat(TextureSemanticClass semanticClass,
                       uint32_t sourceChannelCount,
                       ExrNumericType sourceNumericType,
                       TextureFormat *outFormat,
                       std::string *outError)
{
    if (!outFormat)
    {
        if (outError)
        {
            *outError = "missing pixel format output";
        }
        return false;
    }

    const bool isHalf = (sourceNumericType == ExrNumericType::Float16);
    const bool isFloat = (sourceNumericType == ExrNumericType::Float32);
    if (!isHalf && !isFloat)
    {
        if (outError)
        {
            *outError = "unsupported EXR numeric type for tile format selection";
        }
        return false;
    }

    switch (semanticClass)
    {
        case TextureSemanticClass::Scalar:
            if (sourceChannelCount != 1u)
            {
                if (outError)
                {
                    *outError = "scalar texture must have exactly 1 source channel";
                }
                return false;
            }
            *outFormat =
                isHalf ? TextureFormat::R16_FLOAT : TextureFormat::R32_FLOAT;
            return true;
        case TextureSemanticClass::Normal:
            if (sourceChannelCount < 2u)
            {
                if (outError)
                {
                    *outError = "normal texture must have at least 2 source channels";
                }
                return false;
            }
            *outFormat =
                isHalf ? TextureFormat::RG16_FLOAT : TextureFormat::RG32_FLOAT;
            return true;
        case TextureSemanticClass::Color:
            if (sourceChannelCount != 3u && sourceChannelCount != 4u)
            {
                if (outError)
                {
                    *outError = "color texture must have 3 or 4 source channels";
                }
                return false;
            }
            *outFormat =
                isHalf ? TextureFormat::RGBA16_FLOAT : TextureFormat::RGBA32_FLOAT;
            return true;
        case TextureSemanticClass::Unknown:
        default:
            if (outError)
            {
                *outError = "unknown texture semantic class";
            }
            return false;
    }
}

bool ConvertPixelsToPixelFormat(TextureFormat pixelFormat,
                                const std::vector<float> &rgbaPixels,
                                std::vector<float> *outPixels,
                                std::string *outError)
{
    if (!outPixels)
    {
        if (outError)
        {
            *outError = "missing converted pixel output";
        }
        return false;
    }
    if (!ValidateRgbaPixels(rgbaPixels, outError))
    {
        return false;
    }

    const uint32_t channelCount = ChannelCountForFormat(pixelFormat);
    if (channelCount == 0u)
    {
        if (outError)
        {
            *outError = "invalid virtual texture pixel format";
        }
        return false;
    }

    const size_t pixelCount = rgbaPixels.size() / 4u;
    outPixels->assign(pixelCount * static_cast<size_t>(channelCount), 0.0f);
    for (size_t pixelIndex = 0; pixelIndex < pixelCount; ++pixelIndex)
    {
        const float *src = rgbaPixels.data() + pixelIndex * 4u;
        float *dst = outPixels->data() + pixelIndex * static_cast<size_t>(channelCount);
        switch (channelCount)
        {
            case 1u:
                dst[0] = src[0];
                break;
            case 2u:
                dst[0] = src[0];
                dst[1] = src[1];
                break;
            case 4u:
                std::memcpy(dst, src, 4u * sizeof(float));
                break;
            default:
                if (outError)
                {
                    *outError = "unsupported channel count for conversion";
                }
                return false;
        }
    }
    return true;
}

bool ConvertMipChainToPixelFormat(TextureFormat pixelFormat,
                                  std::vector<tilebin::UdimMipImage> *mipLevels,
                                  std::string *outError)
{
    if (!mipLevels)
    {
        if (outError)
        {
            *outError = "missing mip chain for conversion";
        }
        return false;
    }

    for (tilebin::UdimMipImage &mip : *mipLevels)
    {
        std::vector<float> converted;
        if (!ConvertPixelsToPixelFormat(pixelFormat, mip.rgba, &converted, outError))
        {
            return false;
        }
        mip.pixelFormat = pixelFormat;
        mip.rgba = std::move(converted);
    }
    return true;
}

void ExpandPixelsToRgba(TextureFormat pixelFormat,
                        const std::vector<float> &pixels,
                        std::vector<float> *outRgba)
{
    if (!outRgba)
    {
        return;
    }

    const uint32_t channelCount = ChannelCountForFormat(pixelFormat);
    if (channelCount == 0u || (pixels.size() % static_cast<size_t>(channelCount)) != 0u)
    {
        outRgba->clear();
        return;
    }

    const size_t pixelCount = pixels.size() / static_cast<size_t>(channelCount);
    outRgba->assign(pixelCount * 4u, 0.0f);
    for (size_t pixelIndex = 0; pixelIndex < pixelCount; ++pixelIndex)
    {
        const float *src = pixels.data() + pixelIndex * static_cast<size_t>(channelCount);
        float *dst = outRgba->data() + pixelIndex * 4u;
        switch (channelCount)
        {
            case 1u:
                dst[0] = src[0];
                dst[1] = src[0];
                dst[2] = src[0];
                dst[3] = 1.0f;
                break;
            case 2u:
                dst[0] = src[0];
                dst[1] = src[1];
                dst[2] = 0.0f;
                dst[3] = 1.0f;
                break;
            case 4u:
                std::memcpy(dst, src, 4u * sizeof(float));
                break;
            default:
                outRgba->clear();
                return;
        }
    }
}

bool EncodePixelPayload(TextureFormat pixelFormat,
                        const std::vector<float> &pixels,
                        std::vector<unsigned char> *outBytes,
                        std::string *outError)
{
    if (!outBytes)
    {
        if (outError)
        {
            *outError = "missing encoded payload output";
        }
        return false;
    }

    const uint32_t channelCount = ChannelCountForFormat(pixelFormat);
    const uint32_t bytesPerChannel = TextureFormatBytesPerChannel(pixelFormat);
    if (channelCount == 0u || bytesPerChannel == 0u)
    {
        if (outError)
        {
            *outError = "invalid virtual texture pixel format";
        }
        return false;
    }
    if ((pixels.size() % static_cast<size_t>(channelCount)) != 0u)
    {
        if (outError)
        {
            *outError = "pixel buffer size does not match format";
        }
        return false;
    }

    outBytes->assign(pixels.size() * static_cast<size_t>(bytesPerChannel), 0u);
    if (bytesPerChannel == 4u)
    {
        std::memcpy(outBytes->data(), pixels.data(), outBytes->size());
        return true;
    }

    for (size_t i = 0; i < pixels.size(); ++i)
    {
        const uint16_t halfBits = util::FloatToHalfBits(pixels[i]);
        std::memcpy(outBytes->data() + i * sizeof(uint16_t), &halfBits, sizeof(uint16_t));
    }
    return true;
}

bool DecodePixelPayload(TextureFormat pixelFormat,
                        const unsigned char *bytes,
                        size_t byteCount,
                        std::vector<float> *outPixels,
                        std::string *outError)
{
    if (!outPixels)
    {
        if (outError)
        {
            *outError = "missing decoded pixel output";
        }
        return false;
    }

    const uint32_t channelCount = ChannelCountForFormat(pixelFormat);
    const uint32_t bytesPerChannel = TextureFormatBytesPerChannel(pixelFormat);
    if (channelCount == 0u || bytesPerChannel == 0u)
    {
        if (outError)
        {
            *outError = "invalid virtual texture pixel format";
        }
        return false;
    }
    if ((byteCount % static_cast<size_t>(bytesPerChannel)) != 0u)
    {
        if (outError)
        {
            *outError = "encoded payload size does not match format";
        }
        return false;
    }

    const size_t sampleCount = byteCount / static_cast<size_t>(bytesPerChannel);
    outPixels->assign(sampleCount, 0.0f);
    if (bytesPerChannel == 4u)
    {
        std::memcpy(outPixels->data(), bytes, byteCount);
        return true;
    }

    for (size_t i = 0; i < sampleCount; ++i)
    {
        uint16_t halfBits = 0u;
        std::memcpy(&halfBits, bytes + i * sizeof(uint16_t), sizeof(uint16_t));
        (*outPixels)[i] = util::HalfBitsToFloat(halfBits);
    }
    return true;
}

} // namespace tileprep
} // namespace ybi
