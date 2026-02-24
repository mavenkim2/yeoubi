#include "shared.h"
#include "tile_binary.h"

#include <algorithm>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <unordered_map>
#include <vector>

#include <tinyexr.h>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "third_party/stb_image_write.h"

namespace
{

struct TextureJob
{
    std::string materialPath;
    std::string inputName;
    std::string originalPath;
    std::string resolvedPath;
};

std::string ResolveUdimTilePath(const std::string &path)
{
    const std::string tokenUpper = "<UDIM>";
    const std::string tokenLower = "<udim>";

    size_t pos = path.find(tokenUpper);
    size_t tokenLen = tokenUpper.size();
    if (pos == std::string::npos)
    {
        pos = path.find(tokenLower);
        tokenLen = tokenLower.size();
    }
    if (pos == std::string::npos)
    {
        return path;
    }

    for (int tile = 1001; tile <= 1199; ++tile)
    {
        std::string candidate = path;
        candidate.replace(pos, tokenLen, std::to_string(tile));
        if (std::filesystem::exists(candidate))
        {
            return candidate;
        }
    }

    std::string fallback = path;
    fallback.replace(pos, tokenLen, "1001");
    return fallback;
}

std::string LowerExt(const std::string &path)
{
    std::string ext = std::filesystem::path(path).extension().string();
    for (char &c : ext)
    {
        c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    }
    return ext;
}

bool LoadExrRgba(const std::string &path, int &width, int &height, std::vector<float> &rgba, std::string &reason)
{
    width = 0;
    height = 0;
    auto AssignRgbaAndFree = [&](float *rawData, int w, int h) {
        const size_t pixelCount = static_cast<size_t>(w) * static_cast<size_t>(h);
        rgba.assign(rawData, rawData + pixelCount * 4);
        std::free(rawData);
        width = w;
        height = h;
    };

    float *raw = nullptr;
    const char *err = nullptr;
    int rc = LoadEXR(&raw, &width, &height, path.c_str(), &err);
    if (rc == TINYEXR_SUCCESS && raw)
    {
        AssignRgbaAndFree(raw, width, height);
        return true;
    }

    std::string primaryError;
    if (err)
    {
        primaryError = err;
        FreeEXRErrorMessage(err);
        err = nullptr;
    }

    if (rc != TINYEXR_ERROR_LAYER_NOT_FOUND)
    {
        reason = "LoadEXR failed for " + path;
        if (!primaryError.empty())
        {
            reason += " (" + primaryError + ")";
        }
        return false;
    }

    const char **layerNames = nullptr;
    int numLayers = 0;
    const char *layersErr = nullptr;
    rc = EXRLayers(path.c_str(), &layerNames, &numLayers, &layersErr);
    if (rc != TINYEXR_SUCCESS || !layerNames || numLayers <= 0)
    {
        reason = "LoadEXR failed for " + path;
        if (!primaryError.empty())
        {
            reason += " (" + primaryError + ")";
        }
        if (layersErr)
        {
            reason += " [EXRLayers: ";
            reason += layersErr;
            reason += "]";
            FreeEXRErrorMessage(layersErr);
        }
        if (layerNames)
        {
            for (int i = 0; i < numLayers; ++i)
            {
                if (layerNames[i])
                {
                    std::free((void *)layerNames[i]);
                }
            }
            std::free((void *)layerNames);
        }
        return false;
    }

    std::string layerError;
    for (int i = 0; i < numLayers; ++i)
    {
        float *layerRaw = nullptr;
        int layerWidth = 0;
        int layerHeight = 0;
        const char *layerErr = nullptr;
        const int layerRc =
            LoadEXRWithLayer(&layerRaw, &layerWidth, &layerHeight, path.c_str(), layerNames[i], &layerErr);
        if (layerRc == TINYEXR_SUCCESS && layerRaw)
        {
            AssignRgbaAndFree(layerRaw, layerWidth, layerHeight);
            for (int j = 0; j < numLayers; ++j)
            {
                if (layerNames[j])
                {
                    std::free((void *)layerNames[j]);
                }
            }
            std::free((void *)layerNames);
            return true;
        }
        if (layerErr)
        {
            layerError = layerErr;
            FreeEXRErrorMessage(layerErr);
        }
    }

    for (int i = 0; i < numLayers; ++i)
    {
        if (layerNames[i])
        {
            std::free((void *)layerNames[i]);
        }
    }
    std::free((void *)layerNames);

    reason = "LoadEXRWithLayer failed for " + path;
    if (!layerError.empty())
    {
        reason += " (" + layerError + ")";
    }
    return false;
}

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

void ConvertTileToPngRgba8(const std::vector<float> &tile, std::vector<unsigned char> &outRgba8)
{
    outRgba8.resize(tile.size());
    for (size_t i = 0; i < tile.size(); ++i)
    {
        const float v = std::min(1.0f, std::max(0.0f, tile[i]));
        outRgba8[i] = static_cast<unsigned char>(v * 255.0f + 0.5f);
    }
}

bool WriteTileImagePng(const fs::path &path, int width, int height, const std::vector<float> &tile)
{
    std::vector<unsigned char> rgba8;
    ConvertTileToPngRgba8(tile, rgba8);
    return stbi_write_png(path.string().c_str(), width, height, 4, rgba8.data(), width * 4) != 0;
}

bool WriteTilePreviewImages(const fs::path &verifyDir,
                            const std::string &baseName,
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
            ExtractTileRgbaF32(
                image, imageWidth, imageHeight, tx, ty, tileSize, tilePixels, tileWidth, tileHeight);
            const fs::path verifyPath =
                verifyDir / (baseName + "_tile_" + std::to_string(tx) + "_" + std::to_string(ty) + ".png");
            if (!WriteTileImagePng(verifyPath, tileWidth, tileHeight, tilePixels))
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
                               int expectedWidth,
                               int expectedHeight,
                               const std::vector<float> &sourceRgba,
                               float eps,
                               std::string &outError)
{
    ybi::tilebin::TileFileHeader header = {};
    std::vector<ybi::tilebin::TileRecord> records;
    std::vector<float> reconstructed;
    if (!ybi::tilebin::ReadTileBinary(tilePath, header, records, reconstructed, &outError))
    {
        return false;
    }

    if (static_cast<int>(header.imageWidth) != expectedWidth || static_cast<int>(header.imageHeight) != expectedHeight)
    {
        outError = "tile verify dimension mismatch for " + tilePath.string();
        return false;
    }

    ybi::tilebin::DiffStats diff;
    if (!ybi::tilebin::DiffImagesExact(sourceRgba, reconstructed, eps, &diff))
    {
        outError = "tile verify mismatch for " + tilePath.string() + " mismatches=" +
                   std::to_string(diff.mismatchCount) + " maxAbs=" + std::to_string(diff.maxAbs);
        return false;
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

    std::unordered_map<std::string, TextureJob> uniqueTextures;
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
            TextureJob job = {};
            job.materialPath = mat.materialPath;
            job.inputName = kv.first;
            job.originalPath = kv.second.texturePath;
            job.resolvedPath = ResolveUdimTilePath(kv.second.texturePath);
            uniqueTextures.emplace(job.resolvedPath, std::move(job));
        }
    }

    std::printf("Tile prep: selected materials=%d unique textures=%zu tileSize=%d verifyCount=%d verifyPass=%s eps=%g\n",
                selectedMaterials,
                uniqueTextures.size(),
                cli.tileSize,
                cli.tileVerifyCount,
                cli.tileVerifyPass ? "on" : "off",
                cli.tileVerifyEps);

    int processed = 0;
    for (const auto &it : uniqueTextures)
    {
        const TextureJob &job = it.second;
        if (LowerExt(job.resolvedPath) != ".exr")
        {
            std::printf("Tile prep: skip non-EXR texture %s\n", job.resolvedPath.c_str());
            continue;
        }
        int width = 0;
        int height = 0;
        std::vector<float> rgba;
        std::string reason;
        if (!LoadExrRgba(job.resolvedPath, width, height, rgba, reason))
        {
            if (outError)
            {
                *outError = reason;
            }
            return false;
        }

        const std::string baseName = Sanitize(job.resolvedPath);
        const fs::path binPath = tileOutDir / (baseName + ".tiles.bin");
        if (!ybi::tilebin::WriteTileBinary(binPath,
                                           width,
                                           height,
                                           cli.tileSize,
                                           rgba,
                                           &reason))
        {
            if (outError)
            {
                *outError = reason;
            }
            return false;
        }

        if (!WriteTilePreviewImages(verifyDir,
                                    baseName,
                                    width,
                                    height,
                                    cli.tileSize,
                                    rgba,
                                    cli.tileVerifyCount,
                                    reason))
        {
            if (outError)
            {
                *outError = reason;
            }
            return false;
        }

        if (cli.tileVerifyPass)
        {
            if (!VerifyRoundTripTileBinary(binPath, width, height, rgba, cli.tileVerifyEps, reason))
            {
                if (outError)
                {
                    *outError = reason;
                }
                return false;
            }
        }

        ++processed;
        std::printf("Tile prep: [%d/%zu] %s -> %s (%dx%d)\n",
                    processed,
                    uniqueTextures.size(),
                    job.resolvedPath.c_str(),
                    binPath.string().c_str(),
                    width,
                    height);
    }

    std::printf("Tile prep: done processed=%d\n", processed);
    return true;
}
