#include "core/arena.h"
#include "core/assert.h"
#include "core/half.h"
#include "core/math.h"
#include "core/string_builder.h"
#include "core/thread_context.h"
#include "core/vec4f.h"
#include "io/usd.h"
#include "render/virtual_texture_file.h"
#include "scene/scene.h"
#include "tiff.h"
#include <unordered_map>

#define STB_IMAGE_IMPLEMENTATION
#define TINYEXR_IMPLEMENTATION
#include <algorithm>
#include <atomic>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <limits>
#include <mutex>
#include <stb_image.h>
#include <tbb/parallel_for.h>
#include <tiffio.h>
#include <tinyexr/tinyexr.h>
#include <vector>

namespace ybi
{

struct Image
{
    void *ptr;
    int width;
    int height;
    int numChannels;
    int numBytesPerChannel;
};

static int GetImageFileChannelCount(const std::string &filename)
{
    std::string ext = std::filesystem::path(filename).extension().string();
    std::transform(ext.begin(), ext.end(), ext.begin(), [](unsigned char c) {
        return static_cast<char>(std::tolower(c));
    });

    if (ext == ".png" || ext == ".jpg" || ext == ".jpeg")
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

static void FlipImageVertically(Image *image)
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

static float ReadChannelAsFloat(const void *data,
                                int pixelIndex,
                                int channel,
                                int numChannels,
                                int numBytesPerChannel,
                                SampleFormat sampleFormat)
{
    if (channel >= numChannels)
    {
        return channel == 3 ? 1.0f : 0.0f;
    }

    const int sampleIndex = pixelIndex * numChannels + channel;
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

static void WriteChannelFromFloat(void *data,
                                  int pixelIndex,
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

    const int sampleIndex = pixelIndex * numChannels + channel;
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

static TextureFormat
GetTextureFormat(uint32_t samplesPerPixel, uint32_t bitsPerSample, SampleFormat sampleFormat)
{
    YBI_ASSERT(samplesPerPixel == 1 || samplesPerPixel == 2 || samplesPerPixel == 4);

    if (sampleFormat == SampleFormat::Unorm && bitsPerSample == 8)
    {
        switch (samplesPerPixel)
        {
            case 1:
                return TextureFormat::R8_UNORM;
            case 2:
                return TextureFormat::RG8_UNORM;
            case 4:
                return TextureFormat::RGBA8_UNORM;
        }
    }
    else if (sampleFormat == SampleFormat::Unorm && bitsPerSample == 16)
    {
        switch (samplesPerPixel)
        {
            case 1:
                return TextureFormat::R16_UNORM;
            case 2:
                return TextureFormat::RG16_UNORM;
            case 4:
                return TextureFormat::RGBA16_UNORM;
        }
    }
    else if (sampleFormat == SampleFormat::Float && bitsPerSample == 16)
    {
        switch (samplesPerPixel)
        {
            case 1:
                return TextureFormat::R16_FLOAT;
            case 2:
                return TextureFormat::RG16_FLOAT;
            case 4:
                return TextureFormat::RGBA16_FLOAT;
        }
    }
    else if (sampleFormat == SampleFormat::Float && bitsPerSample == 32)
    {
        switch (samplesPerPixel)
        {
            case 1:
                return TextureFormat::R32_FLOAT;
            case 2:
                return TextureFormat::RG32_FLOAT;
            case 4:
                return TextureFormat::RGBA32_FLOAT;
        }
    }

    YBI_LOGFATAL("Unsupported texture format: %u channels, %u bits, sample format %u\n",
                 samplesPerPixel,
                 bitsPerSample,
                 static_cast<uint32_t>(sampleFormat));
    return TextureFormat::RGBA8_UNORM;
}

static Vec4f CalculateFallbackColor(const Image &texture, SampleFormat sampleFormat)
{
    double sum[4] = {};
    const int numTexels = texture.width * texture.height;

    for (int pixelIndex = 0; pixelIndex < numTexels; pixelIndex++)
    {
        for (int channel = 0; channel < 4; channel++)
        {
            sum[channel] += ReadChannelAsFloat(texture.ptr,
                                               pixelIndex,
                                               channel,
                                               texture.numChannels,
                                               texture.numBytesPerChannel,
                                               sampleFormat);
        }
    }

    const float invNumTexels = numTexels ? 1.0f / static_cast<float>(numTexels) : 0.0f;
    return Vec4f(static_cast<float>(sum[0]) * invNumTexels,
                 static_cast<float>(sum[1]) * invNumTexels,
                 static_cast<float>(sum[2]) * invNumTexels,
                 static_cast<float>(sum[3]) * invNumTexels);
}

static View<uint8_t> GenerateNextMip(Arena &arena,
                                     const void *src,
                                     int srcWidth,
                                     int srcHeight,
                                     int numChannels,
                                     int numBytesPerChannel,
                                     SampleFormat sampleFormat)
{
    const int dstWidth = std::max(1, srcWidth / 2);
    const int dstHeight = std::max(1, srcHeight / 2);
    const int numBytesPerTexel = numChannels * numBytesPerChannel;

    View<uint8_t> dst =
        arena.PushArray<uint8_t>(static_cast<size_t>(dstWidth) * dstHeight * numBytesPerTexel);

    for (int y = 0; y < dstHeight; y++)
    {
        for (int x = 0; x < dstWidth; x++)
        {
            const int dstPixelIndex = y * dstWidth + x;
            for (int channel = 0; channel < numChannels; channel++)
            {
                float sum = 0.0f;
                for (int dy = 0; dy < 2; dy++)
                {
                    for (int dx = 0; dx < 2; dx++)
                    {
                        const int srcX = std::min(srcWidth - 1, x * 2 + dx);
                        const int srcY = std::min(srcHeight - 1, y * 2 + dy);
                        const int srcPixelIndex = srcY * srcWidth + srcX;
                        sum += ReadChannelAsFloat(src,
                                                  srcPixelIndex,
                                                  channel,
                                                  numChannels,
                                                  numBytesPerChannel,
                                                  sampleFormat);
                    }
                }

                WriteChannelFromFloat(dst.data(),
                                      dstPixelIndex,
                                      channel,
                                      numChannels,
                                      numBytesPerChannel,
                                      sampleFormat,
                                      sum * 0.25f);
            }
        }
    }

    return dst;
}

static void WriteTextureMipPages(FileBuilder &fileBuilder,
                                 uint64_t pageTableOffset,
                                 uint32_t firstPageIndex,
                                 const void *mipData,
                                 int mipWidth,
                                 int mipHeight,
                                 int numChannels,
                                 int numBytesPerChannel,
                                 int pageSize,
                                 int gutterSize)
{
    const int numBytesPerTexel = numChannels * numBytesPerChannel;
    const int paddedPageSize = pageSize + 2 * gutterSize;
    const int pageNumBytes = Sqr(paddedPageSize) * numBytesPerTexel;
    const int rowStride = mipWidth * numBytesPerTexel;
    const int pageStride = paddedPageSize * numBytesPerTexel;
    const int xSteps = (mipWidth + pageSize - 1) / pageSize;
    const int ySteps = (mipHeight + pageSize - 1) / pageSize;
    const uint8_t *mipBytes = reinterpret_cast<const uint8_t *>(mipData);

    for (int y = 0; y < ySteps; y++)
    {
        for (int x = 0; x < xSteps; x++)
        {
            const uint32_t pageIndex = firstPageIndex + static_cast<uint32_t>(y * xSteps + x);

            uint64_t offset = fileBuilder.Allocate(pageNumBytes);
            VirtualTexturePageEntry *pageEntry =
                reinterpret_cast<VirtualTexturePageEntry *>(fileBuilder.GetPtr(pageTableOffset)) +
                pageIndex;
            pageEntry->offset = offset;
            pageEntry->packedSize = pageNumBytes;
            pageEntry->unpackedSize = pageNumBytes;

            uint8_t *dst = fileBuilder.GetPtr(offset);
            for (int row = 0; row < paddedPageSize; row++)
            {
                int sourceY = (pageSize * y + row - gutterSize) % mipHeight;
                if (sourceY < 0)
                    sourceY += mipHeight;

                const uint8_t *srcRow = mipBytes + sourceY * rowStride;
                uint8_t *dstRow = dst + row * pageStride;

                int sourceX = (pageSize * x - gutterSize) % mipWidth;
                if (sourceX < 0)
                    sourceX += mipWidth;

                for (int column = 0; column < paddedPageSize;)
                {
                    const int count = Min(paddedPageSize - column, mipWidth - sourceX);
                    memcpy(dstRow + column * numBytesPerTexel,
                           srcRow + sourceX * numBytesPerTexel,
                           count * numBytesPerTexel);
                    column += count;
                    sourceX = 0;
                }
            }
        }
    }
}

static SampleFormat LoadImageFile(Image *texture, const std::string &filename)
{
    int width = 0;
    int height = 0;
    int channels = 0;
    const bool is16Bit = stbi_is_16_bit(filename.c_str()) != 0;
    void *stbiPixels =
        is16Bit
            ? static_cast<void *>(stbi_load_16(filename.c_str(), &width, &height, &channels, 0))
            : static_cast<void *>(stbi_load(filename.c_str(), &width, &height, &channels, 0));
    YBI_ERROR(stbiPixels, "Failed to load image: %s\n", filename.c_str());

    const int numBytesPerChannel = is16Bit ? 2 : 1;
    const size_t numBytes = static_cast<size_t>(width) * height * channels * numBytesPerChannel;
    void *pixels = malloc(numBytes);
    YBI_ASSERT(pixels);
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
    YBI_ERROR(result == TINYEXR_SUCCESS, "Failed to read EXR version: %s\n", filename.c_str());
    YBI_ERROR(!version.multipart, "Multipart EXR is not supported: %s\n", filename.c_str());

    EXRHeader header;
    InitEXRHeader(&header);
    const char *error = nullptr;
    result = ParseEXRHeaderFromFile(&header, &version, filename.c_str(), &error);
    if (result != TINYEXR_SUCCESS)
    {
        if (error)
        {
            printf("TinyEXR error: %s\n", error);
            FreeEXRErrorMessage(error);
        }
        FreeEXRHeader(&header);
        YBI_LOGFATAL("Failed to read EXR header: %s\n", filename.c_str());
    }

    YBI_ERROR(header.num_channels > 0 && header.num_channels <= 4,
              "EXR must contain between one and four channels: %s\n",
              filename.c_str());

    const int pixelType = header.pixel_types[0];
    for (int channel = 1; channel < header.num_channels; channel++)
    {
        YBI_ERROR(header.pixel_types[channel] == pixelType,
                  "Mixed EXR channel data types are not supported: %s\n",
                  filename.c_str());
    }

    YBI_ERROR(pixelType == TINYEXR_PIXELTYPE_HALF || pixelType == TINYEXR_PIXELTYPE_FLOAT,
              "EXR sample type is not supported: %s\n",
              filename.c_str());

    // Keep half channels as half instead of asking TinyEXR to promote them to float.
    for (int channel = 0; channel < header.num_channels; channel++)
    {
        header.requested_pixel_types[channel] = pixelType;
    }

    EXRImage exrImage;
    InitEXRImage(&exrImage);
    error = nullptr;
    result = LoadEXRImageFromFile(&exrImage, &header, filename.c_str(), &error);

    if (result != TINYEXR_SUCCESS)
    {
        if (error)
        {
            printf("TinyEXR error: %s\n", error);
            FreeEXRErrorMessage(error);
        }
        FreeEXRImage(&exrImage);
        FreeEXRHeader(&header);
        YBI_LOGFATAL("Failed to load EXR: %s\n", filename.c_str());
    }

    const int numBytesPerChannel = pixelType == TINYEXR_PIXELTYPE_HALF ? 2 : 4;
    const size_t numTexels = static_cast<size_t>(exrImage.width) * exrImage.height;
    const size_t numBytes = numTexels * header.num_channels * numBytesPerChannel;
    uint8_t *pixels = static_cast<uint8_t *>(malloc(numBytes));
    YBI_ASSERT(pixels);

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
    YBI_ASSERT(channelOrder.size() == static_cast<size_t>(header.num_channels));

    // TinyEXR exposes one plane per channel. Interleave those planes without
    // converting either their values or their scalar representation.
    for (size_t texelIndex = 0; texelIndex < numTexels; texelIndex++)
    {
        for (int channel = 0; channel < header.num_channels; channel++)
        {
            const int srcChannel = channelOrder[channel];
            memcpy(pixels + (texelIndex * header.num_channels + channel) * numBytesPerChannel,
                   exrImage.images[srcChannel] + texelIndex * numBytesPerChannel,
                   numBytesPerChannel);
        }
    }

    texture->ptr = pixels;
    texture->width = exrImage.width;
    texture->height = exrImage.height;
    texture->numChannels = header.num_channels;
    texture->numBytesPerChannel = numBytesPerChannel;

    FreeEXRImage(&exrImage);
    FreeEXRHeader(&header);
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
    TIFF *tiff = TIFFOpen(filename.c_str(), "r");
    YBI_ASSERT(tiff);

    uint32_t width = 0;
    uint32_t height = 0;
    bool result = TIFFGetField(tiff, TIFFTAG_IMAGEWIDTH, &width);
    YBI_ASSERT(result);
    result = TIFFGetField(tiff, TIFFTAG_IMAGELENGTH, &height);
    YBI_ASSERT(result);

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
            YBI_LOGFATAL("Sample format: %u not supported yet\n", sampleFormat);
        }
    }

    YBI_ERROR(samplesPerPixel >= 1 && samplesPerPixel <= 4,
              "TIFF must contain between one and four channels: %s\n",
              filename.c_str());
    YBI_ERROR((textureSampleFormat == SampleFormat::Unorm &&
               (bitsPerSample == 8 || bitsPerSample == 16)) ||
                  (textureSampleFormat == SampleFormat::Float &&
                   (bitsPerSample == 16 || bitsPerSample == 32)),
              "Unsupported TIFF sample type: %u-bit sample format %u\n",
              bitsPerSample,
              sampleFormat);
    YBI_ERROR(planarConfig == PLANARCONFIG_CONTIG || planarConfig == PLANARCONFIG_SEPARATE,
              "Unsupported TIFF planar configuration: %u\n",
              planarConfig);
    YBI_ASSERT(width > 0 && height > 0);
    YBI_ASSERT(width <= static_cast<uint32_t>(std::numeric_limits<int>::max()));
    YBI_ASSERT(height <= static_cast<uint32_t>(std::numeric_limits<int>::max()));

    const int numBytesPerChannel = bitsPerSample / 8;
    const size_t numBytesPerRow =
        static_cast<size_t>(width) * samplesPerPixel * numBytesPerChannel;
    const size_t numBytes = numBytesPerRow * height;
    uint8_t *pixels = static_cast<uint8_t *>(malloc(numBytes));
    YBI_ASSERT(pixels);

    if (TIFFIsTiled(tiff))
    {
        uint32_t tileWidth = 0;
        uint32_t tileHeight = 0;
        YBI_ASSERT(TIFFGetField(tiff, TIFFTAG_TILEWIDTH, &tileWidth));
        YBI_ASSERT(TIFFGetField(tiff, TIFFTAG_TILELENGTH, &tileHeight));

        const tmsize_t tileSize = TIFFTileSize(tiff);
        const tmsize_t tileRowSize = TIFFTileRowSize(tiff);
        YBI_ASSERT(tileSize > 0 && tileRowSize > 0);
        std::vector<uint8_t> tile(static_cast<size_t>(tileSize));

        const uint16_t numPlanes = planarConfig == PLANARCONFIG_CONTIG ? 1 : samplesPerPixel;
        const uint16_t samplesPerPlane = planarConfig == PLANARCONFIG_CONTIG ? samplesPerPixel : 1;
        const size_t sourceTexelSize = samplesPerPlane * numBytesPerChannel;
        YBI_ASSERT(static_cast<size_t>(tileRowSize) >= tileWidth * sourceTexelSize);

        for (uint16_t plane = 0; plane < numPlanes; plane++)
        {
            for (uint32_t tileY = 0; tileY < height; tileY += tileHeight)
            {
                for (uint32_t tileX = 0; tileX < width; tileX += tileWidth)
                {
                    const ttile_t tileIndex = TIFFComputeTile(tiff, tileX, tileY, 0, plane);
                    YBI_ERROR(TIFFReadEncodedTile(tiff, tileIndex, tile.data(), tileSize) >= 0,
                              "Failed to read TIFF tile %u: %s\n",
                              static_cast<uint32_t>(tileIndex),
                              filename.c_str());

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
        YBI_ASSERT(scanlineSize > 0);
        std::vector<uint8_t> scanline(static_cast<size_t>(scanlineSize));

        for (uint32_t row = 0; row < height; row++)
        {
            const uint32_t dstRow = GetTopDownRow(row, height, orientation);
            uint8_t *dst = pixels + static_cast<size_t>(dstRow) * numBytesPerRow;

            if (planarConfig == PLANARCONFIG_CONTIG)
            {
                YBI_ASSERT(static_cast<size_t>(scanlineSize) >= numBytesPerRow);
                YBI_ERROR(TIFFReadScanline(tiff, scanline.data(), row, 0) >= 0,
                          "Failed to read TIFF row %u: %s\n",
                          row,
                          filename.c_str());
                memcpy(dst, scanline.data(), numBytesPerRow);
            }
            else
            {
                const size_t planeRowSize = static_cast<size_t>(width) * numBytesPerChannel;
                YBI_ASSERT(static_cast<size_t>(scanlineSize) >= planeRowSize);
                for (uint16_t channel = 0; channel < samplesPerPixel; channel++)
                {
                    YBI_ERROR(TIFFReadScanline(tiff, scanline.data(), row, channel) >= 0,
                              "Failed to read TIFF row %u, channel %u: %s\n",
                              row,
                              channel,
                              filename.c_str());
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
    TIFFClose(tiff);

    texture->ptr = pixels;
    texture->width = static_cast<int>(width);
    texture->height = static_cast<int>(height);
    texture->numChannels = samplesPerPixel;
    texture->numBytesPerChannel = numBytesPerChannel;
    return textureSampleFormat;
}

static void PadRGBToRGBA(Image *texture, SampleFormat sampleFormat)
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
        WriteChannelFromFloat(rgba,
                              static_cast<int>(texelIndex),
                              3,
                              4,
                              texture->numBytesPerChannel,
                              sampleFormat,
                              1.0f);
    }

    free(texture->ptr);
    texture->ptr = rgba;
    texture->numChannels = 4;
}

static float SrgbToLinear(float value)
{
    return value <= 0.04045f ? value / 12.92f : std::pow((value + 0.055f) / 1.055f, 2.4f);
}

static void LinearizeSrgbTexture(Image *texture, SampleFormat sourceFormat)
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

static void ProcessTexture(Arena &arena,
                           const std::string &filename,
                           uint32_t udim,
                           MaterialColorSpace colorSpace)
{
    std::string ext = std::filesystem::path(filename).extension().string();
    std::transform(ext.begin(), ext.end(), ext.begin(), [](unsigned char c) {
        return static_cast<char>(std::tolower(c));
    });
    Image texture = {};
    SampleFormat sampleFormat = SampleFormat::Count;

    if (ext == ".png" || ext == ".jpg" || ext == ".jpeg")
    {
        sampleFormat = LoadImageFile(&texture, filename);
    }
    else if (ext == ".tif" || ext == ".tiff")
    {
        sampleFormat = LoadTIFFFile(&texture, filename);
    }
    else if (ext == ".exr")
    {
        sampleFormat = LoadEXRFile(&texture, filename);
    }
    else
    {
        YBI_LOGFATAL("File format is not supported: %s\n", filename.c_str());
    }

    PadRGBToRGBA(&texture, sampleFormat);
    if (colorSpace == MaterialColorSpace::SRGB)
    {
        LinearizeSrgbTexture(&texture, sampleFormat);
        sampleFormat = SampleFormat::Float;
    }
    const TextureFormat format =
        GetTextureFormat(texture.numChannels, texture.numBytesPerChannel * 8, sampleFormat);

    // Texture sampling treats v = 0 as the bottom row, so store ybvt source data bottom-up.
    FlipImageVertically(&texture);

    const int pageSize = 128;

    uint32_t numMips = 0;
    uint32_t numPages = 0;
    for (int mipWidth = texture.width, mipHeight = texture.height;
         mipWidth >= pageSize && mipHeight >= pageSize;
         mipWidth = std::max(1, mipWidth / 2), mipHeight = std::max(1, mipHeight / 2))
    {
        const int xSteps = (mipWidth + pageSize - 1) / pageSize;
        const int ySteps = (mipHeight + pageSize - 1) / pageSize;
        numPages += static_cast<uint32_t>(xSteps * ySteps);
        numMips++;
    }

    std::filesystem::path outputPath(filename);
    outputPath.replace_extension(".ybvt");

    FileBuilder fileBuilder(outputPath.string());

    const uint64_t headerOffset = fileBuilder.Allocate(sizeof(VirtualTextureFileHeader));
    const uint64_t pageTableOffset =
        fileBuilder.Allocate(sizeof(VirtualTexturePageEntry) * numPages);

    VirtualTextureFileHeader *header =
        reinterpret_cast<VirtualTextureFileHeader *>(fileBuilder.GetPtr(headerOffset));
    memset(fileBuilder.GetPtr(pageTableOffset), 0, sizeof(VirtualTexturePageEntry) * numPages);

    // TODO: select this based on the filter
    const int gutterSize = 1;

    *header = {};
    header->magic = YBVT_MAGIC;
    header->version = 2;
    header->udim = udim;
    header->pageSize = pageSize;
    header->width = texture.width;
    header->height = texture.height;
    header->numMips = numMips;
    header->gutterSize = gutterSize;
    // header->numChannels = texture.numChannels;
    // header->numBytesPerChannel = texture.numBytesPerChannel;
    header->format = format;
    header->compression = 0;
    header->numPages = numPages;
    header->fallbackColor = CalculateFallbackColor(texture, sampleFormat);
    header->pageTableOffset = pageTableOffset;
    header->pageDataOffset = pageTableOffset + sizeof(VirtualTexturePageEntry) * numPages;

    const void *mipData = texture.ptr;
    int mipWidth = texture.width;
    int mipHeight = texture.height;
    uint32_t pageIndex = 0;

    for (uint32_t mip = 0; mip < numMips; mip++)
    {
        const int xSteps = (mipWidth + pageSize - 1) / pageSize;
        const int ySteps = (mipHeight + pageSize - 1) / pageSize;

        WriteTextureMipPages(fileBuilder,
                             pageTableOffset,
                             pageIndex,
                             mipData,
                             mipWidth,
                             mipHeight,
                             texture.numChannels,
                             texture.numBytesPerChannel,
                             pageSize,
                             gutterSize);
        pageIndex += static_cast<uint32_t>(xSteps * ySteps);

        const int nextMipWidth = std::max(1, mipWidth / 2);
        const int nextMipHeight = std::max(1, mipHeight / 2);
        if (mip + 1 < numMips)
        {
            View<uint8_t> nextMip = GenerateNextMip(arena,
                                                    mipData,
                                                    mipWidth,
                                                    mipHeight,
                                                    texture.numChannels,
                                                    texture.numBytesPerChannel,
                                                    sampleFormat);
            mipData = nextMip.data();
            mipWidth = nextMipWidth;
            mipHeight = nextMipHeight;
        }
    }

    free(texture.ptr);
}

static void Run(const std::string &filename)
{
    DataContainer container;
    Arena arena;
    std::unique_ptr<SceneLoader> sceneLoader = SceneLoader::CreateSceneLoader(filename);
    YBI_ASSERT(sceneLoader);
    SceneImportContext importContext;
    Scene *scene = sceneLoader->LoadScene(&arena, &container, &importContext);
    (void)scene;

    struct TextureToProcess
    {
        std::string filename;
        uint32_t udim;
        MaterialColorSpace colorSpace;
    };
    std::vector<TextureToProcess> textures;
    std::unordered_map<std::string, std::pair<uint32_t, MaterialColorSpace>> textureSettings;
    const auto AddTexture =
        [&](const std::string &filename, uint32_t udim, MaterialColorSpace colorSpace) {
            YBI_ASSERT(!filename.empty());
            const auto [it, inserted] =
                textureSettings.emplace(filename, std::make_pair(udim, colorSpace));
            if (inserted)
            {
                textures.push_back({filename, udim, colorSpace});
            }
            else
            {
                YBI_ERROR(it->second.first == udim,
                          "Texture has conflicting UDIM values: %s (%u and %u)\n",
                          filename.c_str(),
                          it->second.first,
                          udim);
                YBI_ERROR(it->second.second == colorSpace,
                          "Texture is used with conflicting color spaces: %s\n",
                          filename.c_str());
            }
        };

    for (const Material &material : container.materials)
    {
        for (const auto &[diffuseUdimPath, udim] : material.diffuseUdimPaths)
        {
            AddTexture(diffuseUdimPath, udim, MaterialColorSpace::Linear);
        }
    }

#if 0
    std::vector<const UdimTextureVector *> udimTextureVectors;

    for (const OpenPbrMaterial &material : container.openPbrMaterials)
    {
        std::vector<const UdimTextureVector *> materialTextures;
        for (const MaterialNode &node : material.graph.nodes)
        {
            for (const MaterialInput &input : node.inputs)
            {
                const UdimTextureVector *textureFiles =
                    std::get_if<UdimTextureVector>(&input.value);

                if (textureFiles)
                {
                    materialTextures.push_back(textureFiles);
                }
            }
        }

        struct TextureChannel
        {
            std::string filename;
            uint32_t numChannels;
        };

        // Calculate the number of channels per udim
        std::unordered_map<uint32_t, std::vector<TextureChannel>> udimToTextures;
        for (const UdimTextureVector *textureVector : udimTextureVectors)
        {
            // NOTE: assume every texture in a single udim set has the same number of channels
            auto &[file, udim] = (*textureVector)[0];
            uint32_t numChannels = GetImageFileChannelCount(file);

            for (const auto &[filename, udim] : *textureVector)
            {
                auto [it, inserted] = udimToTextures.try_emplace(udim);
                it->second.push_back({filename, numChannels});
            }
        }

        // Iterating over every udim, combine textures into groups of 4 channels if they
        // have the same data type (float, unorm, etc)
        for (const auto &[udim, textures] : udimToTextures)
        {
            size_t numChosen = 0;
            std::vector<bool> chosen;
            std::vector<std::string> combinedTextures;

            while (numChosen != textures.size())
            {
                uint32_t numChannels = 0;
                for (size_t index = 0; index < textures.size(); index++)
                {
                    if (numChannels == 4)
                        break;

                    if (chosen[index])
                        continue;

                    const TextureChannel &channel = textures[index];

                    // TODO IMPORTANT: they have to use the same UV set too
                    uint32_t requestedNumChannels = 4 - numChannels;
                    if (numChannels == 0 || (channel.numChannels == requestedNumChannels))
                    {
                        combinedTextures.push_back(channel.filename);
                        numChannels += channel.numChannels;
                        chosen[index] = true;
                        numChosen++;
                    }
                }

                // TODO: combine textures, somehow have those properly mapped in the renderer
                for (const std::string &textureFile : combinedTextures)
                {
                }

                combinedTextures.clear();
            }
        }

        udimTextureVectors.insert(
            udimTextureVectors.end(), materialTextures.begin(), materialTextures.end());
    }
    for (const UdimTextureVector *textureVector : udimTextureVectors)
    {
        for (const auto &[filename, udim] : *textureVector)
        {
            AddTexture(filename, udim, MaterialColorSpace::Linear);
        }
    }
#endif

    for (const OpenPbrMaterial &material : container.openPbrMaterials)
    {
        for (const MaterialNode &node : material.graph.nodes)
        {
            for (const MaterialInput &input : node.inputs)
            {
                const UdimTextureVector *textureFiles =
                    std::get_if<UdimTextureVector>(&input.value);

                if (textureFiles)
                {
                    for (const auto &[filename, udim] : *textureFiles)
                    {
                        AddTexture(filename, udim, MaterialColorSpace::Linear);
                    }
                }
            }
        }
    }

    for (const MaterialGraph &graph : container.materialXMaterials)
    {
        for (const MaterialNode &node : graph.nodes)
        {
            for (const MaterialInput &input : node.inputs)
            {
                const UdimTextureVector *textureFiles =
                    std::get_if<UdimTextureVector>(&input.value);

                if (textureFiles)
                {
                    for (const auto &[filename, udim] : *textureFiles)
                    {
                        AddTexture(filename, udim, input.colorSpace);
                    }
                }
            }
        }
    }

    const size_t totalTextures = textures.size();
    printf("Num materials: %zu\n", container.materials.size());
    printf("Virtual textures: 0/%zu\n", totalTextures);
    if (totalTextures == 0)
    {
        return;
    }

    std::atomic<size_t> numFinished = 0;
    std::mutex progressMutex;

    tbb::parallel_for(size_t(0), textures.size(), [&](size_t i) {
        ScratchArena scratch;
        const TextureToProcess &texture = textures[i];
        ProcessTexture(*scratch.arena, texture.filename, texture.udim, texture.colorSpace);

        const size_t finished = numFinished.fetch_add(1) + 1;
        {
            std::lock_guard<std::mutex> lock(progressMutex);
            printf("\rVirtual textures: %zu/%zu", finished, totalTextures);
            fflush(stdout);
        }
    });
    printf("\n");
}

} // namespace ybi

using namespace ybi;
int main(int argc, char **argv)
{
    YBI_ASSERT(argc >= 2);
    Run(argv[1]);
    return 0;
}
