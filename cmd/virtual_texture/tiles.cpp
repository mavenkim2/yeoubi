#include "../usd_ntc_encode/shared.h"
#include "tile_binary.h"
#include "texture/exr_io.h"
#include "texture/path_utils.h"
#include "texture/udim_utils.h"

#include <algorithm>
#include <cctype>
#include <cstdio>
#include <cstring>
#include <string>
#include <unordered_map>
#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "third_party/stb_image_write.h"

namespace
{

struct TextureGroup
{
    std::string basePathNoUdim;
    std::unordered_map<uint32_t, std::string> udimPaths;
};

void ExtractTileRgbaF32(const std::vector<float> &image,
                        int imageWidth,
                        int imageHeight,
                        int tileX,
                        int tileY,
                        int tileSize,
                        std::vector<float> &outTile,
                        int &outWidth,
                        int &outHeight)
{
    const int x0 = tileX * tileSize;
    const int y0 = tileY * tileSize;
    outWidth = std::max(0, std::min(tileSize, imageWidth - x0));
    outHeight = std::max(0, std::min(tileSize, imageHeight - y0));
    outTile.assign(static_cast<size_t>(outWidth) * static_cast<size_t>(outHeight) * 4u, 0.0f);
    for (int y = 0; y < outHeight; ++y)
    {
        const float *src = image.data() + (static_cast<size_t>(y0 + y) * static_cast<size_t>(imageWidth) +
                                           static_cast<size_t>(x0)) *
                                              4u;
        float *dst = outTile.data() + static_cast<size_t>(y) * static_cast<size_t>(outWidth) * 4u;
        std::memcpy(dst, src, static_cast<size_t>(outWidth) * 4u * sizeof(float));
    }
}

bool WriteTilePreviewImages(const fs::path &verifyDir,
                            const std::string &baseName,
                            uint32_t udim,
                            int imageWidth,
                            int imageHeight,
                            int tileSize,
                            const std::vector<float> &image,
                            int verifyCount,
                            std::string &outError)
{
    int verifyWritten = 0;
    std::vector<float> tilePixels;
    const int tilesX = (imageWidth + tileSize - 1) / tileSize;
    const int tilesY = (imageHeight + tileSize - 1) / tileSize;
    for (int ty = 0; ty < tilesY && verifyWritten < verifyCount; ++ty)
    {
        for (int tx = 0; tx < tilesX && verifyWritten < verifyCount; ++tx)
        {
            int tileWidth = 0;
            int tileHeight = 0;
            ExtractTileRgbaF32(image, imageWidth, imageHeight, tx, ty, tileSize, tilePixels, tileWidth, tileHeight);
            std::vector<unsigned char> rgba8(tilePixels.size());
            for (size_t i = 0; i < tilePixels.size(); ++i)
            {
                const float v = std::min(1.0f, std::max(0.0f, tilePixels[i]));
                rgba8[i] = static_cast<unsigned char>(v * 255.0f + 0.5f);
            }
            const fs::path verifyPath = verifyDir /
                                        (baseName + "_udim" + std::to_string(udim) + "_tile_" + std::to_string(tx) + "_" +
                                         std::to_string(ty) + ".png");
            if (stbi_write_png(verifyPath.string().c_str(), tileWidth, tileHeight, 4, rgba8.data(), tileWidth * 4) == 0)
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
        ybi::tilebin::DiffStats diff;
        if (!ybi::tilebin::DiffImagesExact(source.rgba, decoded.rgba, eps, &diff))
        {
            outError = "tile verify mismatch udim=" + std::to_string(decoded.udim) +
                       " mismatches=" + std::to_string(diff.mismatchCount) +
                       " maxAbs=" + std::to_string(diff.maxAbs);
            return false;
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
            const std::string &texturePath = kv.second.texturePath;
            const std::string basePathNoUdim = ybi::usd_ntc::StripUdimFromPath(texturePath);
            TextureGroup &group = groups[basePathNoUdim];
            group.basePathNoUdim = basePathNoUdim;
            std::unordered_map<uint32_t, std::string> discovered;
            std::string reason;
            if (!ybi::usd_ntc::CollectUdimPaths(texturePath, discovered, reason))
            {
                if (outError)
                {
                    *outError = reason;
                }
                return false;
            }
            for (auto &entry : discovered)
            {
                group.udimPaths.emplace(entry.first, std::move(entry.second));
            }
        }
    }

    std::printf("Tile prep: selected materials=%d textures=%zu tileSize=%d verifyCount=%d verifyPass=%s eps=%g\n",
                selectedMaterials,
                groups.size(),
                cli.tileSize,
                cli.tileVerifyCount,
                cli.tileVerifyPass ? "on" : "off",
                cli.tileVerifyEps);

    int processed = 0;
    for (const auto &entry : groups)
    {
        const TextureGroup &group = entry.second;
        if (ybi::texture::LowerExt(group.basePathNoUdim) != ".exr")
        {
            std::printf("Tile prep: skip non-EXR texture %s\n", group.basePathNoUdim.c_str());
            continue;
        }

        std::vector<ybi::tilebin::UdimImage> images;
        images.reserve(group.udimPaths.size());
        for (const auto &udimPath : group.udimPaths)
        {
            ybi::tilebin::UdimImage image = {};
            image.udim = udimPath.first;
            int width = 0;
            int height = 0;
            std::string reason;
            if (!ybi::texture::LoadExrRgba(udimPath.second, &width, &height, &image.rgba, &reason, true))
            {
                if (outError)
                {
                    *outError = reason;
                }
                return false;
            }
            image.width = static_cast<uint32_t>(width);
            image.height = static_cast<uint32_t>(height);
            images.push_back(std::move(image));
        }

        const std::string baseName = Sanitize(group.basePathNoUdim);
        const fs::path binPath = tileOutDir / (baseName + ".tiles.bin");
        std::string reason;
        if (!ybi::tilebin::WriteTileBinary(binPath, cli.tileSize, images, &reason))
        {
            if (outError)
            {
                *outError = reason;
            }
            return false;
        }

        for (const ybi::tilebin::UdimImage &img : images)
        {
            if (!WriteTilePreviewImages(verifyDir,
                                        baseName,
                                        img.udim,
                                        static_cast<int>(img.width),
                                        static_cast<int>(img.height),
                                        cli.tileSize,
                                        img.rgba,
                                        cli.tileVerifyCount,
                                        reason))
            {
                if (outError)
                {
                    *outError = reason;
                }
                return false;
            }
        }

        if (cli.tileVerifyPass && !VerifyRoundTripTileBinary(binPath, images, cli.tileVerifyEps, reason))
        {
            if (outError)
            {
                *outError = reason;
            }
            return false;
        }

        ++processed;
        std::printf("Tile prep: [%d/%zu] %s -> %s (udims=%zu)\n",
                    processed,
                    groups.size(),
                    group.basePathNoUdim.c_str(),
                    binPath.string().c_str(),
                    images.size());
    }

    std::printf("Tile prep: done processed=%d\n", processed);
    return true;
}
