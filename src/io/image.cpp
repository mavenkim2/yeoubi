#include "io/image.h"

#include "core/half.h"
#include "core/thread_context.h"
#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <limits>
#include <memory>
#include <utility>

#define STB_IMAGE_IMPLEMENTATION
#include <stb_image.h>
#define TINYEXR_IMPLEMENTATION
#include <tiffio.h>
#include <tinyexr/tinyexr.h>

namespace ybi
{

int GetImageFileChannelCount(const std::string &filename)
{
    std::string ext = std::filesystem::path(filename).extension().string();
    std::transform(ext.begin(), ext.end(), ext.begin(), [](unsigned char c) {
        return static_cast<char>(std::tolower(c));
    });

    if (ext != ".exr" && ext != ".tif" && ext != ".tiff")
    {
        int width = 0;
        int height = 0;
        int numChannels = 0;
        if (!stbi_info(filename.c_str(), &width, &height, &numChannels))
        {
            YBI_LOGFATAL("Failed to read image information: %s\n", filename.c_str());
            return 0;
        }

        return numChannels;
    }

    if (ext == ".tif" || ext == ".tiff")
    {
        TIFF *tiff = TIFFOpen(filename.c_str(), "r");
        if (!tiff)
        {
            YBI_LOGFATAL("Failed to open TIFF: %s\n", filename.c_str());
            return 0;
        }

        uint16_t numChannels = 0;
        const int result = TIFFGetFieldDefaulted(tiff, TIFFTAG_SAMPLESPERPIXEL, &numChannels);
        TIFFClose(tiff);

        if (!result || numChannels == 0)
        {
            YBI_LOGFATAL("Failed to read TIFF channel count: %s\n", filename.c_str());
            return 0;
        }

        return static_cast<int>(numChannels);
    }

    if (ext == ".exr")
    {
        EXRVersion version = {};
        if (ParseEXRVersionFromFile(&version, filename.c_str()) != TINYEXR_SUCCESS)
        {
            YBI_LOGFATAL("Failed to read EXR version: %s\n", filename.c_str());
            return 0;
        }

        if (version.multipart)
        {
            YBI_LOGFATAL("Multipart EXR files do not have a single channel count: %s\n",
                         filename.c_str());
            return 0;
        }

        EXRHeader header;
        InitEXRHeader(&header);

        const char *error = nullptr;
        const int result = ParseEXRHeaderFromFile(&header, &version, filename.c_str(), &error);
        if (result != TINYEXR_SUCCESS)
        {
            if (error)
            {
                printf("TinyEXR error: %s\n", error);
                FreeEXRErrorMessage(error);
            }
            FreeEXRHeader(&header);
            YBI_LOGFATAL("Failed to read EXR header: %s\n", filename.c_str());
            return 0;
        }

        const int numChannels = header.num_channels;
        FreeEXRHeader(&header);

        if (numChannels <= 0)
        {
            YBI_LOGFATAL("EXR has no channels: %s\n", filename.c_str());
            return 0;
        }

        return numChannels;
    }

    YBI_LOGFATAL("Unsupported image file format: %s\n", filename.c_str());
    return 0;
}

void FlipImageVertically(Image *image)
{
    YBI_ASSERT(image);

    if (!image->ptr || image->height <= 1)
    {
        return;
    }

    const size_t rowSize = static_cast<size_t>(image->width) *
                           static_cast<size_t>(image->numChannels) *
                           static_cast<size_t>(image->numBytesPerChannel);
    uint8_t *pixels = reinterpret_cast<uint8_t *>(image->ptr);

    ScratchArena scratchArena;
    View<uint8_t> scratch = scratchArena.arena->PushArray<uint8_t>(rowSize);

    for (int y = 0; y < image->height / 2; y++)
    {
        uint8_t *topRow = pixels + rowSize * static_cast<size_t>(y);
        uint8_t *bottomRow = pixels + rowSize * static_cast<size_t>(image->height - 1 - y);

        memcpy(scratch.data(), topRow, rowSize);
        memcpy(topRow, bottomRow, rowSize);
        memcpy(bottomRow, scratch.data(), rowSize);
    }
}

float ReadChannelAsFloat(const void *data,
                         size_t pixelIndex,
                         int channel,
                         int numChannels,
                         int numBytesPerChannel,
                         SampleFormat sampleFormat)
{
    if (channel >= numChannels)
    {
        return channel == 3 ? 1.0f : 0.0f;
    }

    const size_t sampleIndex = pixelIndex * numChannels + channel;
    if (sampleFormat == SampleFormat::Unorm && numBytesPerChannel == 1)
    {
        const uint8_t *pixels = reinterpret_cast<const uint8_t *>(data);
        return pixels[sampleIndex] / 255.0f;
    }
    else if (sampleFormat == SampleFormat::Unorm && numBytesPerChannel == 2)
    {
        const uint16_t *pixels = reinterpret_cast<const uint16_t *>(data);
        return pixels[sampleIndex] / 65535.0f;
    }
    else if (sampleFormat == SampleFormat::Float && numBytesPerChannel == 2)
    {
        const uint16_t *pixels = reinterpret_cast<const uint16_t *>(data);
        return static_cast<float>(Half::FromBits(pixels[sampleIndex]));
    }
    else if (sampleFormat == SampleFormat::Float && numBytesPerChannel == 4)
    {
        const float *pixels = reinterpret_cast<const float *>(data);
        return pixels[sampleIndex];
    }
    else
    {
        YBI_LOGFATAL("Unsupported sample format and size: %u, %d bytes\n",
                     static_cast<uint32_t>(sampleFormat),
                     numBytesPerChannel);
        return 0.0f;
    }
}

void WriteChannelFromFloat(void *data,
                           size_t pixelIndex,
                           int channel,
                           int numChannels,
                           int numBytesPerChannel,
                           SampleFormat sampleFormat,
                           float value)
{
    if (channel >= numChannels)
    {
        return;
    }

    const size_t sampleIndex = pixelIndex * numChannels + channel;
    if (sampleFormat == SampleFormat::Unorm && numBytesPerChannel == 1)
    {
        uint8_t *pixels = reinterpret_cast<uint8_t *>(data);
        value = std::clamp(value, 0.0f, 1.0f);
        pixels[sampleIndex] = static_cast<uint8_t>(std::round(value * 255.0f));
    }
    else if (sampleFormat == SampleFormat::Unorm && numBytesPerChannel == 2)
    {
        uint16_t *pixels = reinterpret_cast<uint16_t *>(data);
        value = std::clamp(value, 0.0f, 1.0f);
        pixels[sampleIndex] = static_cast<uint16_t>(std::round(value * 65535.0f));
    }
    else if (sampleFormat == SampleFormat::Float && numBytesPerChannel == 2)
    {
        uint16_t *pixels = reinterpret_cast<uint16_t *>(data);
        pixels[sampleIndex] = Half(value).Bits();
    }
    else if (sampleFormat == SampleFormat::Float && numBytesPerChannel == 4)
    {
        float *pixels = reinterpret_cast<float *>(data);
        pixels[sampleIndex] = value;
    }
    else
    {
        YBI_LOGFATAL("Unsupported sample format and size: %u, %d bytes\n",
                     static_cast<uint32_t>(sampleFormat),
                     numBytesPerChannel);
    }
}

static SampleFormat LoadStbImageFile(Image *texture, const std::string &filename)
{
    int width = 0;
    int height = 0;
    int channels = 0;
    if (stbi_is_hdr(filename.c_str()))
    {
        float *stbiPixels = stbi_loadf(filename.c_str(), &width, &height, &channels, 0);
        if (!stbiPixels)
            return SampleFormat::Count;

        const size_t numBytes = size_t(width) * height * channels * sizeof(float);
        void *pixels = malloc(numBytes);
        if (!pixels)
        {
            stbi_image_free(stbiPixels);
            return SampleFormat::Count;
        }
        memcpy(pixels, stbiPixels, numBytes);
        stbi_image_free(stbiPixels);

        texture->ptr = pixels;
        texture->width = width;
        texture->height = height;
        texture->numChannels = channels;
        texture->numBytesPerChannel = sizeof(float);
        return SampleFormat::Float;
    }

    const bool is16Bit = stbi_is_16_bit(filename.c_str()) != 0;
    void *stbiPixels =
        is16Bit
            ? static_cast<void *>(stbi_load_16(filename.c_str(), &width, &height, &channels, 0))
            : static_cast<void *>(stbi_load(filename.c_str(), &width, &height, &channels, 0));
    if (!stbiPixels)
        return SampleFormat::Count;

    const int numBytesPerChannel = is16Bit ? 2 : 1;
    const size_t numBytes = static_cast<size_t>(width) * height * channels * numBytesPerChannel;
    void *pixels = malloc(numBytes);
    if (!pixels)
    {
        stbi_image_free(stbiPixels);
        return SampleFormat::Count;
    }
    memcpy(pixels, stbiPixels, numBytes);
    stbi_image_free(stbiPixels);

    texture->ptr = pixels;
    texture->width = width;
    texture->height = height;
    texture->numChannels = channels;
    texture->numBytesPerChannel = numBytesPerChannel;
    return SampleFormat::Unorm;
}

static bool EXRChannelHasComponent(const char *name, char component)
{
    const size_t length = strlen(name);
    return length > 0 && name[length - 1] == component && (length == 1 || name[length - 2] == '.');
}

static SampleFormat LoadEXRFile(Image *texture, const std::string &filename)
{
    EXRVersion version = {};
    int result = ParseEXRVersionFromFile(&version, filename.c_str());
    if (result != TINYEXR_SUCCESS || version.multipart)
        return SampleFormat::Count;

    EXRHeader header;
    InitEXRHeader(&header);
    std::unique_ptr<EXRHeader, decltype(&FreeEXRHeader)> headerStorage(&header, &FreeEXRHeader);
    const char *error = nullptr;
    result = ParseEXRHeaderFromFile(&header, &version, filename.c_str(), &error);
    if (result != TINYEXR_SUCCESS)
    {
        if (error)
        {
            printf("TinyEXR error: %s\n", error);
            FreeEXRErrorMessage(error);
        }
        return SampleFormat::Count;
    }

    if (header.num_channels <= 0 || header.num_channels > 4)
        return SampleFormat::Count;

    const int pixelType = header.pixel_types[0];
    for (int channel = 1; channel < header.num_channels; channel++)
    {
        if (header.pixel_types[channel] != pixelType)
            return SampleFormat::Count;
    }

    if (pixelType != TINYEXR_PIXELTYPE_HALF && pixelType != TINYEXR_PIXELTYPE_FLOAT)
        return SampleFormat::Count;

    // Keep half channels as half instead of asking TinyEXR to promote them to float.
    for (int channel = 0; channel < header.num_channels; channel++)
    {
        header.requested_pixel_types[channel] = pixelType;
    }

    EXRImage exrImage;
    InitEXRImage(&exrImage);
    std::unique_ptr<EXRImage, decltype(&FreeEXRImage)> imageStorage(&exrImage, &FreeEXRImage);
    error = nullptr;
    result = LoadEXRImageFromFile(&exrImage, &header, filename.c_str(), &error);

    if (result != TINYEXR_SUCCESS)
    {
        if (error)
        {
            printf("TinyEXR error: %s\n", error);
            FreeEXRErrorMessage(error);
        }
        return SampleFormat::Count;
    }

    if (exrImage.width <= 0 || exrImage.height <= 0 || (!exrImage.images && !exrImage.tiles))
        return SampleFormat::Count;

    const int numBytesPerChannel = pixelType == TINYEXR_PIXELTYPE_HALF ? 2 : 4;
    const size_t numTexels = static_cast<size_t>(exrImage.width) * exrImage.height;
    const size_t numBytes = numTexels * header.num_channels * numBytesPerChannel;
    uint8_t *pixels = static_cast<uint8_t *>(malloc(numBytes));
    std::unique_ptr<uint8_t, decltype(&std::free)> pixelStorage(pixels, &std::free);
    if (!pixels)
        return SampleFormat::Count;

    std::vector<int> channelOrder;
    const auto AppendEXRComponent = [&](char component) {
        for (int channel = 0; channel < header.num_channels; channel++)
        {
            if (EXRChannelHasComponent(header.channels[channel].name, component))
            {
                channelOrder.push_back(channel);
                return;
            }
        }
    };

    // TinyEXR commonly exposes standard color channels in A, B, G, R order.
    // Put recognizable color or luminance channels into texture component order.
    int numColorComponents = 0;
    int numLuminanceComponents = 0;
    for (int channel = 0; channel < header.num_channels; channel++)
    {
        const char *name = header.channels[channel].name;
        numColorComponents +=
            EXRChannelHasComponent(name, 'R') || EXRChannelHasComponent(name, 'G') ||
                    EXRChannelHasComponent(name, 'B') || EXRChannelHasComponent(name, 'A')
                ? 1
                : 0;
        numLuminanceComponents +=
            EXRChannelHasComponent(name, 'Y') || EXRChannelHasComponent(name, 'A') ? 1 : 0;
    }

    if (numColorComponents == header.num_channels)
    {
        AppendEXRComponent('R');
        AppendEXRComponent('G');
        AppendEXRComponent('B');
        AppendEXRComponent('A');
    }
    else if (numLuminanceComponents == header.num_channels)
    {
        AppendEXRComponent('Y');
        AppendEXRComponent('A');
    }
    else
    {
        for (int channel = 0; channel < header.num_channels; channel++)
        {
            channelOrder.push_back(channel);
        }
    }
    if (channelOrder.size() != static_cast<size_t>(header.num_channels))
        return SampleFormat::Count;

    // TinyEXR exposes one plane per channel. Interleave those planes without
    // converting either their values or their scalar representation.
    const auto CopyTexel = [&](unsigned char *const *channels, size_t srcIndex, size_t dstIndex) {
        for (int channel = 0; channel < header.num_channels; channel++)
        {
            const int srcChannel = channelOrder[channel];
            memcpy(pixels + (dstIndex * header.num_channels + channel) * numBytesPerChannel,
                   channels[srcChannel] + srcIndex * numBytesPerChannel,
                   numBytesPerChannel);
        }
    };

    if (exrImage.images)
    {
        for (size_t texelIndex = 0; texelIndex < numTexels; texelIndex++)
            CopyTexel(exrImage.images, texelIndex, texelIndex);
    }
    else
    {
        if (header.tile_size_x <= 0 || header.tile_size_y <= 0 || exrImage.num_tiles <= 0)
            return SampleFormat::Count;
        for (int tileIndex = 0; tileIndex < exrImage.num_tiles; tileIndex++)
        {
            const EXRTile &tile = exrImage.tiles[tileIndex];
            const size_t startX = size_t(tile.offset_x) * header.tile_size_x;
            const size_t startY = size_t(tile.offset_y) * header.tile_size_y;
            if (!tile.images || tile.offset_x < 0 || tile.offset_y < 0 || tile.width <= 0 ||
                tile.height <= 0 || tile.width > header.tile_size_x ||
                tile.height > header.tile_size_y || startX + tile.width > size_t(exrImage.width) ||
                startY + tile.height > size_t(exrImage.height))
                return SampleFormat::Count;
            for (int y = 0; y < tile.height; y++)
            {
                for (int x = 0; x < tile.width; x++)
                {
                    const size_t srcIndex = size_t(y) * header.tile_size_x + x;
                    const size_t dstIndex = (startY + y) * exrImage.width + startX + x;
                    CopyTexel(tile.images, srcIndex, dstIndex);
                }
            }
        }
    }

    texture->ptr = pixelStorage.release();
    texture->width = exrImage.width;
    texture->height = exrImage.height;
    texture->numChannels = header.num_channels;
    texture->numBytesPerChannel = numBytesPerChannel;

    return SampleFormat::Float;
}

static uint32_t GetTopDownRow(uint32_t row, uint32_t height, uint16_t orientation)
{
    if (orientation == ORIENTATION_TOPLEFT)
    {
        return row;
    }
    if (orientation == ORIENTATION_BOTLEFT)
    {
        return height - 1u - row;
    }

    YBI_LOGFATAL("Unsupported TIFF orientation: %u\n", orientation);
    return 0;
}

static SampleFormat LoadTIFFFile(Image *texture, const std::string &filename)
{
    std::unique_ptr<TIFF, decltype(&TIFFClose)> file(TIFFOpen(filename.c_str(), "r"), &TIFFClose);
    TIFF *tiff = file.get();
    if (!tiff)
        return SampleFormat::Count;

    uint32_t width = 0;
    uint32_t height = 0;
    if (!TIFFGetField(tiff, TIFFTAG_IMAGEWIDTH, &width) ||
        !TIFFGetField(tiff, TIFFTAG_IMAGELENGTH, &height))
        return SampleFormat::Count;

    uint16_t samplesPerPixel = 0;
    uint16_t bitsPerSample = 0;
    uint16_t sampleFormat = 0;
    uint16_t planarConfig = 0;
    uint16_t orientation = 0;

    TIFFGetFieldDefaulted(tiff, TIFFTAG_SAMPLESPERPIXEL, &samplesPerPixel);
    TIFFGetFieldDefaulted(tiff, TIFFTAG_BITSPERSAMPLE, &bitsPerSample);
    TIFFGetFieldDefaulted(tiff, TIFFTAG_SAMPLEFORMAT, &sampleFormat);
    TIFFGetFieldDefaulted(tiff, TIFFTAG_PLANARCONFIG, &planarConfig);
    TIFFGetFieldDefaulted(tiff, TIFFTAG_ORIENTATION, &orientation);

    SampleFormat textureSampleFormat;
    switch (sampleFormat)
    {
        case SAMPLEFORMAT_UINT:
        {
            // TODO: this assumes that the uints are interpreted as unorms.
            // need a more robust way of identifying unorm vs actual uints
            textureSampleFormat = SampleFormat::Unorm;
            break;
        }
        case SAMPLEFORMAT_IEEEFP:
        {
            textureSampleFormat = SampleFormat::Float;
            break;
        }
        default:
        {
            return SampleFormat::Count;
        }
    }

    const bool supportedSampleType = (textureSampleFormat == SampleFormat::Unorm &&
                                      (bitsPerSample == 8 || bitsPerSample == 16)) ||
                                     (textureSampleFormat == SampleFormat::Float &&
                                      (bitsPerSample == 16 || bitsPerSample == 32));
    if (samplesPerPixel < 1 || samplesPerPixel > 4 || !supportedSampleType ||
        (planarConfig != PLANARCONFIG_CONTIG && planarConfig != PLANARCONFIG_SEPARATE) ||
        (orientation != ORIENTATION_TOPLEFT && orientation != ORIENTATION_BOTLEFT) || width == 0 ||
        height == 0 || width > static_cast<uint32_t>(std::numeric_limits<int>::max()) ||
        height > static_cast<uint32_t>(std::numeric_limits<int>::max()))
        return SampleFormat::Count;

    const int numBytesPerChannel = bitsPerSample / 8;
    const size_t numBytesPerRow =
        static_cast<size_t>(width) * samplesPerPixel * numBytesPerChannel;
    const size_t numBytes = numBytesPerRow * height;
    uint8_t *pixels = static_cast<uint8_t *>(malloc(numBytes));
    std::unique_ptr<uint8_t, decltype(&std::free)> pixelStorage(pixels, &std::free);
    if (!pixels)
        return SampleFormat::Count;

    if (TIFFIsTiled(tiff))
    {
        uint32_t tileWidth = 0;
        uint32_t tileHeight = 0;
        if (!TIFFGetField(tiff, TIFFTAG_TILEWIDTH, &tileWidth) ||
            !TIFFGetField(tiff, TIFFTAG_TILELENGTH, &tileHeight) || tileWidth == 0 ||
            tileHeight == 0)
            return SampleFormat::Count;

        const tmsize_t tileSize = TIFFTileSize(tiff);
        const tmsize_t tileRowSize = TIFFTileRowSize(tiff);
        if (tileSize <= 0 || tileRowSize <= 0)
            return SampleFormat::Count;
        std::vector<uint8_t> tile(static_cast<size_t>(tileSize));

        const uint16_t numPlanes = planarConfig == PLANARCONFIG_CONTIG ? 1 : samplesPerPixel;
        const uint16_t samplesPerPlane = planarConfig == PLANARCONFIG_CONTIG ? samplesPerPixel : 1;
        const size_t sourceTexelSize = samplesPerPlane * numBytesPerChannel;
        if (static_cast<size_t>(tileRowSize) < tileWidth * sourceTexelSize)
            return SampleFormat::Count;

        for (uint16_t plane = 0; plane < numPlanes; plane++)
        {
            for (uint32_t tileY = 0; tileY < height; tileY += tileHeight)
            {
                for (uint32_t tileX = 0; tileX < width; tileX += tileWidth)
                {
                    const ttile_t tileIndex = TIFFComputeTile(tiff, tileX, tileY, 0, plane);
                    if (TIFFReadEncodedTile(tiff, tileIndex, tile.data(), tileSize) < 0)
                        return SampleFormat::Count;

                    const uint32_t copyWidth = std::min(tileWidth, width - tileX);
                    const uint32_t copyHeight = std::min(tileHeight, height - tileY);
                    for (uint32_t localY = 0; localY < copyHeight; localY++)
                    {
                        const uint32_t dstRow = GetTopDownRow(tileY + localY, height, orientation);
                        uint8_t *dstRowBytes =
                            pixels + static_cast<size_t>(dstRow) * numBytesPerRow;
                        const uint8_t *srcRowBytes =
                            tile.data() + static_cast<size_t>(localY) * tileRowSize;

                        for (uint32_t localX = 0; localX < copyWidth; localX++)
                        {
                            uint8_t *dst = dstRowBytes + static_cast<size_t>(tileX + localX) *
                                                             samplesPerPixel * numBytesPerChannel;
                            const uint8_t *src =
                                srcRowBytes + static_cast<size_t>(localX) * sourceTexelSize;
                            if (planarConfig == PLANARCONFIG_CONTIG)
                            {
                                memcpy(dst, src, sourceTexelSize);
                            }
                            else
                            {
                                memcpy(dst + plane * numBytesPerChannel, src, numBytesPerChannel);
                            }
                        }
                    }
                }
            }
        }
    }
    else
    {
        const tmsize_t scanlineSize = TIFFScanlineSize(tiff);
        if (scanlineSize <= 0)
            return SampleFormat::Count;
        std::vector<uint8_t> scanline(static_cast<size_t>(scanlineSize));

        for (uint32_t row = 0; row < height; row++)
        {
            const uint32_t dstRow = GetTopDownRow(row, height, orientation);
            uint8_t *dst = pixels + static_cast<size_t>(dstRow) * numBytesPerRow;

            if (planarConfig == PLANARCONFIG_CONTIG)
            {
                if (static_cast<size_t>(scanlineSize) < numBytesPerRow ||
                    TIFFReadScanline(tiff, scanline.data(), row, 0) < 0)
                    return SampleFormat::Count;
                memcpy(dst, scanline.data(), numBytesPerRow);
            }
            else
            {
                const size_t planeRowSize = static_cast<size_t>(width) * numBytesPerChannel;
                if (static_cast<size_t>(scanlineSize) < planeRowSize)
                    return SampleFormat::Count;
                for (uint16_t channel = 0; channel < samplesPerPixel; channel++)
                {
                    if (TIFFReadScanline(tiff, scanline.data(), row, channel) < 0)
                        return SampleFormat::Count;
                    for (uint32_t x = 0; x < width; x++)
                    {
                        memcpy(dst + (static_cast<size_t>(x) * samplesPerPixel + channel) *
                                         numBytesPerChannel,
                               scanline.data() + static_cast<size_t>(x) * numBytesPerChannel,
                               numBytesPerChannel);
                    }
                }
            }
        }
    }
    texture->ptr = pixelStorage.release();
    texture->width = static_cast<int>(width);
    texture->height = static_cast<int>(height);
    texture->numChannels = samplesPerPixel;
    texture->numBytesPerChannel = numBytesPerChannel;
    return textureSampleFormat;
}

void PadRGBToRGBA(Image *texture, SampleFormat sampleFormat)
{
    if (texture->numChannels != 3)
    {
        return;
    }

    const size_t numTexels = static_cast<size_t>(texture->width) * texture->height;
    const size_t srcTexelSize = 3u * texture->numBytesPerChannel;
    const size_t dstTexelSize = 4u * texture->numBytesPerChannel;
    void *rgba = malloc(numTexels * dstTexelSize);
    YBI_ASSERT(rgba);

    const uint8_t *src = static_cast<const uint8_t *>(texture->ptr);
    uint8_t *dst = static_cast<uint8_t *>(rgba);
    for (size_t texelIndex = 0; texelIndex < numTexels; texelIndex++)
    {
        memcpy(dst + texelIndex * dstTexelSize, src + texelIndex * srcTexelSize, srcTexelSize);
        WriteChannelFromFloat(
            rgba, texelIndex, 3, 4, texture->numBytesPerChannel, sampleFormat, 1.0f);
    }

    free(texture->ptr);
    texture->ptr = rgba;
    texture->numChannels = 4;
}

static float SrgbToLinear(float value)
{
    return value <= 0.04045f ? value / 12.92f : std::pow((value + 0.055f) / 1.055f, 2.4f);
}

void LinearizeSrgbTexture(Image *texture, SampleFormat sourceFormat)
{
    const int numTexels = texture->width * texture->height;
    const int colorChannels = texture->numChannels == 2 ? 1 : std::min(texture->numChannels, 3);
    const int linearBytesPerChannel =
        sourceFormat == SampleFormat::Float && texture->numBytesPerChannel == 4 ? 4 : 2;
    void *linearPixels =
        malloc(static_cast<size_t>(numTexels) * texture->numChannels * linearBytesPerChannel);
    YBI_ASSERT(linearPixels);

    for (int pixelIndex = 0; pixelIndex < numTexels; pixelIndex++)
    {
        for (int channel = 0; channel < texture->numChannels; channel++)
        {
            float value = ReadChannelAsFloat(texture->ptr,
                                             pixelIndex,
                                             channel,
                                             texture->numChannels,
                                             texture->numBytesPerChannel,
                                             sourceFormat);
            if (channel < colorChannels)
                value = SrgbToLinear(value);

            WriteChannelFromFloat(linearPixels,
                                  pixelIndex,
                                  channel,
                                  texture->numChannels,
                                  linearBytesPerChannel,
                                  SampleFormat::Float,
                                  value);
        }
    }

    free(texture->ptr);
    texture->ptr = linearPixels;
    texture->numBytesPerChannel = linearBytesPerChannel;
}

SampleFormat LoadImageFile(Image *image, const std::string &filename)
{
    std::string extension = std::filesystem::path(filename).extension().string();
    std::transform(extension.begin(), extension.end(), extension.begin(), [](unsigned char c) {
        return static_cast<char>(std::tolower(c));
    });

    SampleFormat format;
    if (extension == ".exr")
        format = LoadEXRFile(image, filename);
    else if (extension == ".tif" || extension == ".tiff")
        format = LoadTIFFFile(image, filename);
    else
        format = LoadStbImageFile(image, filename);

    if (format == SampleFormat::Count)
        std::fprintf(stderr, "Failed to load image: %s\n", filename.c_str());
    return format;
}

bool LoadLinearRGBImage(const std::string &filename,
                        uint32_t &width,
                        uint32_t &height,
                        std::vector<Vec3f> &pixels)
{
    Image image;
    const SampleFormat format = LoadImageFile(&image, filename);
    std::unique_ptr<void, decltype(&std::free)> storage(image.ptr, &std::free);
    if (format == SampleFormat::Count)
        return false;

    const size_t numPixels = size_t(image.width) * image.height;
    std::vector<Vec3f> rgb(numPixels);
    for (size_t pixelIndex = 0; pixelIndex < numPixels; pixelIndex++)
    {
        const float r = ReadChannelAsFloat(
            image.ptr, pixelIndex, 0, image.numChannels, image.numBytesPerChannel, format);
        const float g = image.numChannels <= 2 ? r
                                               : ReadChannelAsFloat(image.ptr,
                                                                    pixelIndex,
                                                                    1,
                                                                    image.numChannels,
                                                                    image.numBytesPerChannel,
                                                                    format);
        const float b = image.numChannels <= 2 ? r
                                               : ReadChannelAsFloat(image.ptr,
                                                                    pixelIndex,
                                                                    2,
                                                                    image.numChannels,
                                                                    image.numBytesPerChannel,
                                                                    format);
        rgb[pixelIndex] = format == SampleFormat::Unorm
                              ? Vec3f(SrgbToLinear(r), SrgbToLinear(g), SrgbToLinear(b))
                              : Vec3f(r, g, b);
    }

    width = static_cast<uint32_t>(image.width);
    height = static_cast<uint32_t>(image.height);
    pixels = std::move(rgb);
    return true;
}

} // namespace ybi
