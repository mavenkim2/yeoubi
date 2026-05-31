#include "shared.h"
#include "exr_mips.h"
#include "tile_binary_detail.h"
#include "texture/path_utils.h"
#include "texture/udim_utils.h"
#include "tile_binary.h"
#include "tile_pixels.h"

#include <algorithm>
#include <atomic>
#include <cctype>
#include <cstdio>
#include <cstring>
#include <functional>
#include <mutex>
#include <string>
#include <unordered_map>
#include <vector>

#include <tbb/parallel_for.h>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "third_party/stb_image_write.h"

namespace
{

struct TextureGroup
{
    std::string basePathNoUdim;
    std::string texturePath;
    std::string firstInputName;
    std::string conflictingInputName;
    ybi::tileprep::TextureSemanticClass semanticClass = ybi::tileprep::TextureSemanticClass::Unknown;
};

ybi::tileprep::TextureSemanticClass ClassifyInputSemantic(const std::string &inputName)
{
    if (inputName == "roughness" || inputName == "metallic" || inputName == "occlusion" ||
        inputName == "ior" || inputName == "opacity" || inputName == "clearcoat" ||
        inputName == "clearcoatRoughness")
    {
        return ybi::tileprep::TextureSemanticClass::Scalar;
    }
    if (inputName == "normal")
    {
        return ybi::tileprep::TextureSemanticClass::Normal;
    }
    return ybi::tileprep::TextureSemanticClass::Color;
}

bool AccumulateTextureSemantic(TextureGroup *group, const std::string &inputName)
{
    if (!group)
    {
        return false;
    }
    const ybi::tileprep::TextureSemanticClass semanticClass = ClassifyInputSemantic(inputName);
    if (group->semanticClass == ybi::tileprep::TextureSemanticClass::Unknown)
    {
        group->semanticClass = semanticClass;
        group->firstInputName = inputName;
        return true;
    }
    if (group->semanticClass == semanticClass)
    {
        return true;
    }
    if (group->conflictingInputName.empty())
    {
        group->conflictingInputName = inputName;
    }
    return false;
}

std::string MakeTileOutputStem(const std::string &basePathNoUdim)
{
    return std::to_string(
               static_cast<unsigned long long>(std::hash<std::string>{}(basePathNoUdim))) +
           "_" + Sanitize(fs::path(basePathNoUdim).filename().string());
}

bool WriteTilePreviewImages(const fs::path &verifyDir,
                            const std::string &baseName,
                            uint32_t udim,
                            int imageWidth,
                            int imageHeight,
                            int tileSize,
                            ybi::TextureFormat pixelFormat,
                            const std::vector<float> &image,
                            int verifyCount,
                            std::string &outError)
{
    int verifyWritten = 0;
    std::vector<float> rgbaImage;
    ybi::tileprep::ExpandPixelsToRgba(pixelFormat, image, &rgbaImage);
    if (rgbaImage.empty())
    {
        outError = "failed expanding preview image to RGBA";
        return false;
    }
    std::vector<float> tilePixels;
    const int tilesX = (imageWidth + tileSize - 1) / tileSize;
    const int tilesY = (imageHeight + tileSize - 1) / tileSize;
    for (int ty = 0; ty < tilesY && verifyWritten < verifyCount; ++ty)
    {
        for (int tx = 0; tx < tilesX && verifyWritten < verifyCount; ++tx)
        {
            int tileWidth = 0;
            int tileHeight = 0;
            ybi::tilebin::detail::ExtractTileFloatSamples(rgbaImage,
                                                          imageWidth,
                                                          imageHeight,
                                                          4u,
                                                          tx,
                                                          ty,
                                                          tileSize,
                                                          tilePixels,
                                                          tileWidth,
                                                          tileHeight);
            std::vector<unsigned char> rgba8(tilePixels.size());
            for (size_t i = 0; i < tilePixels.size(); ++i)
            {
                const float v = std::min(1.0f, std::max(0.0f, tilePixels[i]));
                rgba8[i] = static_cast<unsigned char>(v * 255.0f + 0.5f);
            }
            const fs::path verifyPath =
                verifyDir / (baseName + "_udim" + std::to_string(udim) + "_tile_" +
                             std::to_string(tx) + "_" + std::to_string(ty) + ".png");
            if (stbi_write_png(verifyPath.string().c_str(),
                               tileWidth,
                               tileHeight,
                               4,
                               rgba8.data(),
                               tileWidth * 4) == 0)
            {
                outError = "failed writing tile verify image: " + verifyPath.string();
                return false;
            }
            ++verifyWritten;
        }
    }
    return true;
}

bool VerifyRoundTripTileBinary(const fs::path &tilePath,
                               const std::vector<ybi::tilebin::UdimImage> &sourceImages,
                               float eps,
                               std::string &outError)
{
    ybi::tilebin::TileFileHeader header = {};
    std::vector<ybi::tilebin::UdimEntry> entries;
    std::vector<ybi::tilebin::UdimImage> reconstructed;
    if (!ybi::tilebin::ReadTileBinary(tilePath, header, entries, reconstructed, &outError))
    {
        return false;
    }
    if (sourceImages.size() != reconstructed.size())
    {
        outError = "tile verify UDIM count mismatch for " + tilePath.string();
        return false;
    }

    std::unordered_map<uint32_t, const ybi::tilebin::UdimImage *> sourceByUdim;
    for (const ybi::tilebin::UdimImage &img : sourceImages)
    {
        sourceByUdim.emplace(img.udim, &img);
    }

    for (const ybi::tilebin::UdimImage &decoded : reconstructed)
    {
        const auto it = sourceByUdim.find(decoded.udim);
        if (it == sourceByUdim.end())
        {
            outError = "tile verify missing source udim=" + std::to_string(decoded.udim);
            return false;
        }
            const ybi::tilebin::UdimImage &source = *it->second;
            if (source.width != decoded.width || source.height != decoded.height)
        {
            outError = "tile verify dimension mismatch udim=" + std::to_string(decoded.udim);
            return false;
        }

        if (source.mipLevels.size() != decoded.mipLevels.size())
        {
            outError = "tile verify mip count mismatch udim=" + std::to_string(decoded.udim);
            return false;
        }

        for (size_t mipIndex = 0; mipIndex < source.mipLevels.size(); ++mipIndex)
        {
            const ybi::tilebin::UdimMipImage &sourceMip = source.mipLevels[mipIndex];
            const ybi::tilebin::UdimMipImage &decodedMip = decoded.mipLevels[mipIndex];
            if (sourceMip.level != decodedMip.level || sourceMip.width != decodedMip.width ||
                sourceMip.height != decodedMip.height || sourceMip.pixelFormat != decodedMip.pixelFormat)
            {
                outError =
                    "tile verify mip metadata mismatch udim=" + std::to_string(decoded.udim) +
                    " mip=" + std::to_string(mipIndex);
                return false;
            }

            ybi::tilebin::DiffStats diff;
            if (!ybi::tilebin::DiffImagesExact(sourceMip.rgba, decodedMip.rgba, eps, &diff))
            {
                outError = "tile verify mismatch udim=" + std::to_string(decoded.udim) +
                           " mip=" + std::to_string(mipIndex) +
                           " mismatches=" + std::to_string(diff.mismatchCount) +
                           " maxAbs=" + std::to_string(diff.maxAbs);
                return false;
            }
        }
    }
    return true;
}

} // namespace

bool PrepareTexturesForStreamingTiles(const std::vector<MaterialChannels> &materials,
                                      const Cli &cli,
                                      std::string *outError)
{
    const fs::path tileOutDir = fs::path(cli.outDir) / "stream_tiles";
    const fs::path verifyDir = tileOutDir / "verify_tiles";
    std::error_code ec;
    fs::create_directories(tileOutDir, ec);
    if (ec)
    {
        if (outError)
        {
            *outError = "failed to create stream tile out dir: " + tileOutDir.string();
        }
        return false;
    }
    fs::create_directories(verifyDir, ec);
    if (ec)
    {
        if (outError)
        {
            *outError = "failed to create verify tile out dir: " + verifyDir.string();
        }
        return false;
    }

    std::unordered_map<std::string, TextureGroup> groups;
    int selectedMaterials = 0;
    for (const MaterialChannels &mat : materials)
    {
        if (mat.channels.empty())
        {
            continue;
        }
        if (cli.maxMaterials > 0 && selectedMaterials >= cli.maxMaterials)
        {
            break;
        }
        ++selectedMaterials;

        for (const auto &kv : mat.channels)
        {
            const std::string &inputName = kv.first;
            const std::string &texturePath = kv.second.texturePath;
            const std::string basePathNoUdim = ybi::texture::StripUdimFromPath(texturePath);
            TextureGroup &group = groups[basePathNoUdim];
            group.basePathNoUdim = basePathNoUdim;
            if (group.texturePath.empty())
            {
                group.texturePath = texturePath;
            }
            if (!AccumulateTextureSemantic(&group, inputName))
            {
                if (outError)
                {
                    *outError = "Tile prep: texture path " + basePathNoUdim +
                                " is used by incompatible semantics (" + group.firstInputName +
                                " vs " + group.conflictingInputName + ")";
                }
                return false;
            }
        }
    }

    const int tileSize = 128;
    if (cli.tileSize != tileSize)
    {
        std::printf("Tile prep: forcing tileSize=%d (requested=%d)\n", tileSize, cli.tileSize);
    }

    std::printf("Tile prep: selected materials=%d textures=%zu tileSize=%d verifyCount=%d "
                "verifyPass=%s eps=%g\n",
                selectedMaterials,
                groups.size(),
                tileSize,
                cli.tileVerifyCount,
                cli.tileVerifyPass ? "on" : "off",
                cli.tileVerifyEps);

    std::vector<const TextureGroup *> workItems;
    workItems.reserve(groups.size());
    for (const auto &entry : groups)
    {
        workItems.push_back(&entry.second);
    }

    std::atomic<int> processed{0};
    std::atomic<int> failed{0};
    std::mutex logMutex;

    tbb::parallel_for(size_t(0), workItems.size(), [&](size_t i) {
        const TextureGroup &group = *workItems[i];
        if (ybi::texture::LowerExt(group.basePathNoUdim) != ".exr")
        {
            std::lock_guard<std::mutex> lock(logMutex);
            std::printf("Tile prep: skip non-EXR texture %s\n", group.basePathNoUdim.c_str());
            return;
        }

        std::vector<std::pair<uint32_t, std::string>> udimPaths;
        std::string reason;
        if (!ybi::texture::CollectUdimPaths(group.texturePath, udimPaths, reason))
        {
            if (reason.rfind("missing texture file ", 0) == 0 ||
                reason.rfind("no UDIM files found for ", 0) == 0)
            {
                std::lock_guard<std::mutex> lock(logMutex);
                std::printf("Tile prep: skip missing texture %s\n", group.texturePath.c_str());
                return;
            }
            failed.fetch_add(1, std::memory_order_relaxed);
            std::lock_guard<std::mutex> lock(logMutex);
            std::printf("Tile prep: FAIL %s : %s\n", group.basePathNoUdim.c_str(), reason.c_str());
            return;
        }

        std::vector<ybi::tilebin::UdimImage> images;
        images.reserve(udimPaths.size());
        bool groupFailed = false;
        bool havePixelFormat = false;
        ybi::TextureFormat pixelFormat = ybi::TextureFormat::RGBA32_FLOAT;
        for (const auto &udimPath : udimPaths)
        {
            ybi::tilebin::UdimImage image = {};
            image.udim = udimPath.first;

            ybi::tileprep::ExrMipChainResult mipResult = {};
            std::string mipReason;
            if (!ybi::tileprep::LoadExrMipChain(udimPath.second, true, &mipResult, &mipReason))
            {
                failed.fetch_add(1, std::memory_order_relaxed);
                std::lock_guard<std::mutex> lock(logMutex);
                std::printf(
                    "Tile prep: FAIL %s : %s\n", group.basePathNoUdim.c_str(), mipReason.c_str());
                groupFailed = true;
                break;
            }

            {
                std::lock_guard<std::mutex> lock(logMutex);
                std::printf("Tile prep: exr %s udim=%u mipmap_levels=%s tiled=%d level_mode=%d "
                            "loaded_mips=%zu\n",
                            udimPath.second.c_str(),
                            udimPath.first,
                            mipResult.hasStoredMipLevels ? "yes" : "no",
                            mipResult.tiled,
                            mipResult.tileLevelMode,
                            mipResult.mipLevels.size());
            }

            if (mipResult.mipLevels.empty())
            {
                failed.fetch_add(1, std::memory_order_relaxed);
                std::lock_guard<std::mutex> lock(logMutex);
                std::printf("Tile prep: FAIL %s : empty mip chain for %s\n",
                            group.basePathNoUdim.c_str(),
                            udimPath.second.c_str());
                groupFailed = true;
                break;
            }

            image.width = mipResult.mipLevels[0].width;
            image.height = mipResult.mipLevels[0].height;
            ybi::TextureFormat udimPixelFormat = ybi::TextureFormat::RGBA32_FLOAT;
            if (!ybi::tileprep::ChoosePixelFormat(group.semanticClass,
                                                  mipResult.sourceChannelCount,
                                                  mipResult.sourceNumericType,
                                                  &udimPixelFormat,
                                                  &mipReason))
            {
                failed.fetch_add(1, std::memory_order_relaxed);
                std::lock_guard<std::mutex> lock(logMutex);
                std::printf(
                    "Tile prep: FAIL %s : %s\n", group.basePathNoUdim.c_str(), mipReason.c_str());
                groupFailed = true;
                break;
            }
            if (havePixelFormat && udimPixelFormat != pixelFormat)
            {
                failed.fetch_add(1, std::memory_order_relaxed);
                std::lock_guard<std::mutex> lock(logMutex);
                std::printf("Tile prep: FAIL %s : inconsistent pixel format across UDIMs\n",
                            group.basePathNoUdim.c_str());
                groupFailed = true;
                break;
            }
            pixelFormat = udimPixelFormat;
            havePixelFormat = true;
            if (!ybi::tileprep::ConvertMipChainToPixelFormat(
                    pixelFormat, &mipResult.mipLevels, &mipReason))
            {
                failed.fetch_add(1, std::memory_order_relaxed);
                std::lock_guard<std::mutex> lock(logMutex);
                std::printf(
                    "Tile prep: FAIL %s : %s\n", group.basePathNoUdim.c_str(), mipReason.c_str());
                groupFailed = true;
                break;
            }
            image.pixelFormat = pixelFormat;
            image.mipLevels = std::move(mipResult.mipLevels);
            images.push_back(std::move(image));
        }
        if (groupFailed)
        {
            return;
        }

        const std::string baseName = MakeTileOutputStem(group.basePathNoUdim);
        const fs::path binPath = tileOutDir / (baseName + ".tiles.bin");
        if (!ybi::tilebin::WriteTileBinary(binPath, tileSize, images, &reason))
        {
            failed.fetch_add(1, std::memory_order_relaxed);
            std::lock_guard<std::mutex> lock(logMutex);
            std::printf("Tile prep: FAIL %s : %s\n", group.basePathNoUdim.c_str(), reason.c_str());
            return;
        }

        for (const ybi::tilebin::UdimImage &img : images)
        {
            if (!WriteTilePreviewImages(verifyDir,
                                        baseName,
                                        img.udim,
                                        static_cast<int>(img.width),
                                        static_cast<int>(img.height),
                                        tileSize,
                                        img.pixelFormat,
                                        img.mipLevels[0].rgba,
                                        cli.tileVerifyCount,
                                        reason))
            {
                failed.fetch_add(1, std::memory_order_relaxed);
                std::lock_guard<std::mutex> lock(logMutex);
                std::printf(
                    "Tile prep: FAIL %s : %s\n", group.basePathNoUdim.c_str(), reason.c_str());
                groupFailed = true;
                break;
            }
        }
        if (groupFailed)
        {
            return;
        }

        if (cli.tileVerifyPass &&
            !VerifyRoundTripTileBinary(binPath, images, cli.tileVerifyEps, reason))
        {
            failed.fetch_add(1, std::memory_order_relaxed);
            std::lock_guard<std::mutex> lock(logMutex);
            std::printf("Tile prep: FAIL %s : %s\n", group.basePathNoUdim.c_str(), reason.c_str());
            return;
        }

        const int done = processed.fetch_add(1, std::memory_order_relaxed) + 1;
        {
            std::lock_guard<std::mutex> lock(logMutex);
            std::printf("Tile prep: [%d/%zu] %s -> %s (udims=%zu)\n",
                        done,
                        groups.size(),
                        group.basePathNoUdim.c_str(),
                        binPath.string().c_str(),
                        images.size());
        }
    });

    std::printf("Tile prep: done processed=%d failed=%d total=%zu\n",
                processed.load(std::memory_order_relaxed),
                failed.load(std::memory_order_relaxed),
                groups.size());
    return true;
}
