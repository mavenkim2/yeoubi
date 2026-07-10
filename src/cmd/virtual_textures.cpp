#include "core/arena.h"
#include "core/assert.h"
#include "core/math.h"
#include "core/string_builder.h"
#include "core/thread_context.h"
#include "core/vec4f.h"
#include "io/usd.h"
#include "render/virtual_texture_file.h"
#include "scene/scene.h"
#include "scene/texture.h"
#include <unordered_set>

#define STB_IMAGE_IMPLEMENTATION
#define TINYEXR_IMPLEMENTATION
#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <mutex>
#include <stb_image.h>
#include <tbb/parallel_for.h>
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

static float ReadChannelAsFloat(
    const void *data, int pixelIndex, int channel, int numChannels, int numBytesPerChannel)
{
    if (channel >= numChannels)
    {
        return channel == 3 ? 1.0f : 0.0f;
    }

    if (numBytesPerChannel == 1)
    {
        const uint8_t *pixels = reinterpret_cast<const uint8_t *>(data);
        return pixels[pixelIndex * numChannels + channel] / 255.0f;
    }
    else if (numBytesPerChannel == 4)
    {
        const float *pixels = reinterpret_cast<const float *>(data);
        return pixels[pixelIndex * numChannels + channel];
    }
    else
    {
        YBI_LOGFATAL("Unsupported bytes per channel: %d\n", numBytesPerChannel);
        return 0.0f;
    }
}

static void WriteChannelFromFloat(
    void *data, int pixelIndex, int channel, int numChannels, int numBytesPerChannel, float value)
{
    if (channel >= numChannels)
    {
        return;
    }

    if (numBytesPerChannel == 1)
    {
        uint8_t *pixels = reinterpret_cast<uint8_t *>(data);
        value = std::clamp(value, 0.0f, 1.0f);
        pixels[pixelIndex * numChannels + channel] =
            static_cast<uint8_t>(std::round(value * 255.0f));
    }
    else if (numBytesPerChannel == 4)
    {
        float *pixels = reinterpret_cast<float *>(data);
        pixels[pixelIndex * numChannels + channel] = value;
    }
    else
    {
        YBI_LOGFATAL("Unsupported bytes per channel: %d\n", numBytesPerChannel);
    }
}

static Vec4f CalculateFallbackColor(const Image &texture)
{
    double sum[4] = {};
    const int numTexels = texture.width * texture.height;

    for (int pixelIndex = 0; pixelIndex < numTexels; pixelIndex++)
    {
        for (int channel = 0; channel < 4; channel++)
        {
            sum[channel] += ReadChannelAsFloat(
                texture.ptr, pixelIndex, channel, texture.numChannels, texture.numBytesPerChannel);
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
                                     int numBytesPerChannel)
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
                        sum += ReadChannelAsFloat(
                            src, srcPixelIndex, channel, numChannels, numBytesPerChannel);
                    }
                }

                WriteChannelFromFloat(dst.data(),
                                      dstPixelIndex,
                                      channel,
                                      numChannels,
                                      numBytesPerChannel,
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
                                 int pageSize)
{
    const int numBytesPerTexel = numChannels * numBytesPerChannel;
    const int pageNumBytes = pageSize * pageSize * numBytesPerTexel;
    const int rowStride = mipWidth * numBytesPerTexel;
    const int pageStride = pageSize * numBytesPerTexel;
    const int xSteps = (mipWidth + pageSize - 1) / pageSize;
    const int ySteps = (mipHeight + pageSize - 1) / pageSize;

    for (int y = 0; y < ySteps; y++)
    {
        for (int x = 0; x < xSteps; x++)
        {
            const uint32_t pageIndex = firstPageIndex + static_cast<uint32_t>(y * xSteps + x);

            const int yStart = pageSize * y;
            const int yEnd = Min(mipHeight, pageSize * (y + 1));

            const int xStart = pageSize * x;
            const int xEnd = Min(mipWidth, pageSize * (x + 1));

            const int sizeToCopy = (xEnd - xStart) * numBytesPerTexel;

            uint64_t offset = fileBuilder.Allocate(pageNumBytes);
            VirtualTexturePageEntry *pageEntry =
                reinterpret_cast<VirtualTexturePageEntry *>(fileBuilder.GetPtr(pageTableOffset)) +
                pageIndex;
            pageEntry->offset = offset;
            pageEntry->packedSize = pageNumBytes;
            pageEntry->unpackedSize = pageNumBytes;

            uint8_t *dst = fileBuilder.GetPtr(offset);
            memset(dst, 0, pageNumBytes);

            const uint8_t *src = reinterpret_cast<const uint8_t *>(mipData) +
                                 (yStart * mipWidth + xStart) * numBytesPerTexel;

            for (int row = yStart; row < yEnd; row++)
            {
                memcpy(dst, src, sizeToCopy);
                dst += pageStride;
                src += rowStride;
            }
        }
    }
}

static void LoadImageFile(Image *texture, std::string filename)
{
    int width, height, channels;
    unsigned char *buffer = stbi_load(filename.c_str(), &width, &height, &channels, 4);
    YBI_ASSERT(buffer);

    texture->ptr = reinterpret_cast<void *>(buffer);
    texture->width = width;
    texture->height = height;
    texture->numChannels = 4;
    texture->numBytesPerChannel = 1;
}

static void LoadEXRFile(Image *texture, std::string filename)
{
    float *rgba = nullptr;
    int width = 0;
    int height = 0;
    const char *err = nullptr;

    int result = LoadEXR(&rgba, &width, &height, filename.c_str(), &err);

    if (result != TINYEXR_SUCCESS)
    {
        if (err)
        {
            FreeEXRErrorMessage(err);
            err = nullptr;
        }

        const char **layerNames = nullptr;
        int numLayers = 0;
        result = EXRLayers(filename.c_str(), &layerNames, &numLayers, &err);

        if (result == TINYEXR_SUCCESS && numLayers > 0)
        {
            result =
                LoadEXRWithLayer(&rgba, &width, &height, filename.c_str(), layerNames[0], &err);
        }

        if (layerNames)
        {
            for (int layerIndex = 0; layerIndex < numLayers; layerIndex++)
            {
                free(const_cast<char *>(layerNames[layerIndex]));
            }
            free(layerNames);
        }

        if (result != TINYEXR_SUCCESS)
        {
            if (err)
            {
                printf("TinyEXR error: %s\n", err);
                FreeEXRErrorMessage(err);
            }
            YBI_LOGFATAL("Failed to load EXR: %s\n", filename.c_str());
        }
    }

    texture->ptr = reinterpret_cast<void *>(rgba);
    texture->width = width;
    texture->height = height;
    texture->numChannels = 4;
    texture->numBytesPerChannel = 4;
}

static void ProcessTexture(Arena &arena, const std::string &filename)
{
    std::string ext = std::filesystem::path(filename).extension().string();
    Image texture;
    TextureFormat format;

    if (ext.find("png") != std::string::npos || ext.find("jpg") != std::string::npos)
    {
        LoadImageFile(&texture, filename);
        format = TextureFormat::RGBA8_UNORM;
    }
    else if (ext.find("exr") != std::string::npos)
    {
        LoadEXRFile(&texture, filename);
        format = TextureFormat::RGBA32_FLOAT;
    }
    else
    {
        YBI_LOGFATAL("File format is not supported: %s\n", filename.c_str());
    }

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

    *header = {};
    header->magic = YBVT_MAGIC;
    header->version = 1;
    header->pageSize = pageSize;
    header->width = texture.width;
    header->height = texture.height;
    header->numMips = numMips;
    // header->numChannels = texture.numChannels;
    // header->numBytesPerChannel = texture.numBytesPerChannel;
    header->format = format;
    header->compression = 0;
    header->numPages = numPages;
    header->fallbackColor = CalculateFallbackColor(texture);
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
                             pageSize);
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
                                                    texture.numBytesPerChannel);
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
    Scene *scene = sceneLoader->LoadScene(&arena, &container);
    (void)scene;

    std::vector<std::string> filenames;
    std::unordered_map<std::string, std::string> sss;
    for (const Material &material : container.materials)
    {
        for (const auto &[diffuseUdimPath, udim] : material.diffuseUdimPaths)
        {
            YBI_ASSERT(!diffuseUdimPath.empty());
            auto it = sss.find(diffuseUdimPath);
            if (it == sss.end())
            {
                filenames.push_back(diffuseUdimPath);
                sss.insert({diffuseUdimPath, material.path});
            }
            else
            {
                // printf("%s %s %s\n\n",
                //        diffuseUdimPath.c_str(),
                //        material.path.c_str(),
                //        it->second.c_str());
            }
        }
    }

    const size_t totalTextures = filenames.size();
    printf("Num materials: %zu\n", container.materials.size());
    printf("Virtual textures: 0/%zu\n", totalTextures);
    if (totalTextures == 0)
    {
        return;
    }

    std::atomic<size_t> numFinished = 0;
    std::mutex progressMutex;

    tbb::parallel_for(size_t(0), filenames.size(), [&](size_t i) {
        ScratchArena scratch;
        ProcessTexture(*scratch.arena, filenames[i]);

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
