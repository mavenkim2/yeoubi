#include "core/arena.h"
#include "core/assert.h"
#include "core/math.h"
#include "core/string_builder.h"
#include "core/thread_context.h"
#include "core/vec4f.h"
#include "io/image.h"
#include "io/usd.h"
#include "render/virtual_texture/file.h"
#include "scene/scene.h"
#include <algorithm>
#include <atomic>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <mutex>
#include <tbb/parallel_for.h>
#include <unordered_map>
#include <vector>

namespace ybi
{

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
    const int xSteps = CeilDiv(mipWidth, pageSize);
    const int ySteps = CeilDiv(mipHeight, pageSize);
    const uint8_t *mipBytes = reinterpret_cast<const uint8_t *>(mipData);
    YBI_ASSERT(gutterSize <= Min(mipWidth, mipHeight));

    for (int y = 0; y < ySteps; y++)
    {
        for (int x = 0; x < xSteps; x++)
        {
            const uint32_t pageIndex = firstPageIndex + static_cast<uint32_t>(y * xSteps + x);
            const int yStart = Max(0, pageSize * y - gutterSize);
            const int yEnd = Min(mipHeight, pageSize * (y + 1) + gutterSize);
            const int xStart = Max(0, pageSize * x - gutterSize);
            const int xEnd = Min(mipWidth, pageSize * (x + 1) + gutterSize);
            const int copyWidth = (xEnd - xStart) * numBytesPerTexel;
            const int dstX = x == 0 ? gutterSize : 0;
            const int dstY = y == 0 ? gutterSize : 0;

            uint64_t offset = fileBuilder.Allocate(pageNumBytes);
            VirtualTexturePageEntry *pageEntry =
                reinterpret_cast<VirtualTexturePageEntry *>(fileBuilder.GetPtr(pageTableOffset)) +
                pageIndex;
            pageEntry->offset = offset;
            pageEntry->packedSize = pageNumBytes;
            pageEntry->unpackedSize = pageNumBytes;

            uint8_t *dst = fileBuilder.GetPtr(offset);
            memset(dst, 0, pageNumBytes);

            auto copyRow = [&](int sourceY, int destinationY) {
                const uint8_t *srcRow = mipBytes + sourceY * rowStride;
                uint8_t *dstRow = dst + destinationY * pageStride;
                memcpy(dstRow + dstX * numBytesPerTexel,
                       srcRow + xStart * numBytesPerTexel,
                       copyWidth);

                if (gutterSize && x == 0)
                {
                    memcpy(dstRow,
                           srcRow + (mipWidth - gutterSize) * numBytesPerTexel,
                           gutterSize * numBytesPerTexel);
                }
                if (gutterSize && x == xSteps - 1)
                {
                    memcpy(dstRow + dstX * numBytesPerTexel + copyWidth,
                           srcRow,
                           gutterSize * numBytesPerTexel);
                }
            };

            for (int row = yStart; row < yEnd; row++)
            {
                copyRow(row, dstY + row - yStart);
            }

            if (gutterSize && y == 0)
            {
                for (int row = 0; row < gutterSize; row++)
                    copyRow(mipHeight - gutterSize + row, row);
            }
            if (gutterSize && y == ySteps - 1)
            {
                for (int row = 0; row < gutterSize; row++)
                    copyRow(row, dstY + yEnd - yStart + row);
            }
        }
    }
}

static void ProcessTexture(Arena &arena,
                           const std::string &filename,
                           uint32_t udim,
                           MaterialColorSpace colorSpace)
{
    Image texture;
    SampleFormat sampleFormat = LoadImageFile(&texture, filename);
    if (sampleFormat == SampleFormat::Count)
        return;

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
        const int xSteps = CeilDiv(mipWidth, pageSize);
        const int ySteps = CeilDiv(mipHeight, pageSize);
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
        const int xSteps = CeilDiv(mipWidth, pageSize);
        const int ySteps = CeilDiv(mipHeight, pageSize);

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
