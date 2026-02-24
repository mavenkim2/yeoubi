#include "tile_binary.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <unordered_map>
#include <vector>

#define TINYEXR_IMPLEMENTATION
#include <tinyexr.h>

namespace
{

std::string Sanitize(const std::string &s)
{
    std::string out = s;
    for (char &c : out)
    {
        if (!(std::isalnum((unsigned char)c) || c == '_' || c == '-' || c == '.'))
        {
            c = '_';
        }
    }
    if (out.empty())
    {
        out = "texture";
    }
    return out;
}

bool FindUdimToken(const std::string &path, size_t &pos, size_t &len)
{
    pos = path.find("<UDIM>");
    len = 6;
    if (pos != std::string::npos)
    {
        return true;
    }
    pos = path.find("<udim>");
    len = 6;
    return pos != std::string::npos;
}

bool TryFindUdimDigits(const std::string &path, uint32_t &udim, size_t &digitPos)
{
    const size_t ext = path.rfind('.');
    if (ext == std::string::npos || ext < 4)
    {
        return false;
    }
    digitPos = ext - 4;
    for (size_t i = 0; i < 4; ++i)
    {
        if (!std::isdigit(static_cast<unsigned char>(path[digitPos + i])))
        {
            return false;
        }
    }
    if (digitPos > 0)
    {
        const char c = path[digitPos - 1];
        if (!(c == '.' || c == '_' || c == '-'))
        {
            return false;
        }
    }
    udim = static_cast<uint32_t>(std::strtoul(path.substr(digitPos, 4).c_str(), nullptr, 10));
    return udim >= 1001 && udim <= 1999;
}

std::string StripUdimFromPath(const std::string &path)
{
    size_t tokenPos = std::string::npos;
    size_t tokenLen = 0;
    if (FindUdimToken(path, tokenPos, tokenLen))
    {
        std::string out = path;
        out.erase(tokenPos, tokenLen);
        if (tokenPos > 0 && tokenPos < out.size() && out[tokenPos - 1] == '.' && out[tokenPos] == '.')
        {
            out.erase(tokenPos - 1, 1);
        }
        if (tokenPos > 0 && tokenPos < out.size() &&
            (out[tokenPos - 1] == '_' || out[tokenPos - 1] == '-') && out[tokenPos] == '.')
        {
            out.erase(tokenPos - 1, 1);
        }
        return out;
    }

    uint32_t udim = 0;
    size_t digitPos = 0;
    if (TryFindUdimDigits(path, udim, digitPos))
    {
        std::string out = path;
        out.erase(digitPos, 4);
        if (digitPos > 0 && digitPos < out.size() && out[digitPos - 1] == '.' && out[digitPos] == '.')
        {
            out.erase(digitPos - 1, 1);
        }
        if (digitPos > 0 && digitPos < out.size() &&
            (out[digitPos - 1] == '_' || out[digitPos - 1] == '-') && out[digitPos] == '.')
        {
            out.erase(digitPos - 1, 1);
        }
        return out;
    }

    return path;
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

void PrintUsage(const char *exe)
{
    std::printf("Usage:\n");
    std::printf("  %s --exr <path.exr> --tiles-bin <path.tiles.bin> [--udim <n>] [--eps <float>]\n", exe);
    std::printf("  %s --exr <path.exr> --tiles-dir <dir> [--udim <n>] [--eps <float>]\n", exe);
    std::printf("Notes:\n");
    std::printf("  --tiles-dir auto-resolves to <dir>/<sanitize(strip_udim(exr))>.tiles.bin\n");
}

} // namespace

int main(int argc, char **argv)
{
    std::string exrPath;
    std::string tilesBinPath;
    std::string tilesDir;
    float eps = 1e-7f;
    uint32_t udim = 0;
    bool udimProvided = false;

    for (int i = 1; i < argc; ++i)
    {
        const std::string arg = argv[i];
        auto readValue = [&](std::string &dst) -> bool {
            if (i + 1 >= argc)
            {
                return false;
            }
            dst = argv[++i];
            return true;
        };
        if (arg == "--exr")
        {
            if (!readValue(exrPath))
            {
                PrintUsage(argv[0]);
                return 1;
            }
        }
        else if (arg == "--tiles-bin")
        {
            if (!readValue(tilesBinPath))
            {
                PrintUsage(argv[0]);
                return 1;
            }
        }
        else if (arg == "--tiles-dir")
        {
            if (!readValue(tilesDir))
            {
                PrintUsage(argv[0]);
                return 1;
            }
        }
        else if (arg == "--eps")
        {
            std::string value;
            if (!readValue(value))
            {
                PrintUsage(argv[0]);
                return 1;
            }
            eps = std::max(0.0f, std::strtof(value.c_str(), nullptr));
        }
        else if (arg == "--udim")
        {
            std::string value;
            if (!readValue(value))
            {
                PrintUsage(argv[0]);
                return 1;
            }
            udim = static_cast<uint32_t>(std::strtoul(value.c_str(), nullptr, 10));
            udimProvided = true;
        }
        else if (arg == "--help" || arg == "-h")
        {
            PrintUsage(argv[0]);
            return 0;
        }
        else
        {
            std::printf("Unknown option: %s\n", arg.c_str());
            PrintUsage(argv[0]);
            return 1;
        }
    }

    if (exrPath.empty())
    {
        PrintUsage(argv[0]);
        return 1;
    }

    if (!udimProvided)
    {
        size_t digitPos = 0;
        if (!TryFindUdimDigits(exrPath, udim, digitPos))
        {
            udim = 1001;
        }
    }

    if (tilesBinPath.empty())
    {
        if (tilesDir.empty())
        {
            std::printf("missing tiles input: set --tiles-bin or --tiles-dir\n");
            return 1;
        }
        tilesBinPath = tilesDir + "/" + Sanitize(StripUdimFromPath(exrPath)) + ".tiles.bin";
    }

    int exrW = 0;
    int exrH = 0;
    std::vector<float> exrRgba;
    std::string error;
    if (!LoadExrRgba(exrPath, exrW, exrH, exrRgba, error))
    {
        std::printf("verify failed: %s\n", error.c_str());
        return 2;
    }

    ybi::tilebin::TileFileHeader header = {};
    std::vector<ybi::tilebin::UdimEntry> entries;
    std::vector<ybi::tilebin::UdimImage> images;
    if (!ybi::tilebin::ReadTileBinary(tilesBinPath, header, entries, images, &error))
    {
        std::printf("verify failed: %s\n", error.c_str());
        return 2;
    }

    const ybi::tilebin::UdimImage *decoded = nullptr;
    for (const ybi::tilebin::UdimImage &img : images)
    {
        if (img.udim == udim)
        {
            decoded = &img;
            break;
        }
    }
    if (!decoded)
    {
        std::printf("verify failed: udim %u not found in %s\n", udim, tilesBinPath.c_str());
        return 2;
    }
    if (static_cast<int>(decoded->width) != exrW || static_cast<int>(decoded->height) != exrH)
    {
        std::printf("verify failed: dimension mismatch exr=%dx%d tile=%ux%u udim=%u\n",
                    exrW,
                    exrH,
                    decoded->width,
                    decoded->height,
                    udim);
        return 2;
    }

    ybi::tilebin::DiffStats diff = {};
    const bool ok = ybi::tilebin::DiffImagesExact(exrRgba, decoded->rgba, eps, &diff);
    std::printf("verify stats: exr=%s tiles=%s udim=%u\n", exrPath.c_str(), tilesBinPath.c_str(), udim);
    std::printf("verify stats: width=%d height=%d channels=4\n", exrW, exrH);
    std::printf("verify stats: eps=%g maxAbs=%.9g meanAbs=%.9g rmse=%.9g mismatches=%llu/%zu\n",
                eps,
                diff.maxAbs,
                diff.meanAbs,
                diff.rmse,
                static_cast<unsigned long long>(diff.mismatchCount),
                exrRgba.size());

    if (!ok && diff.firstMismatch < exrRgba.size())
    {
        const size_t pixel = diff.firstMismatch / 4u;
        const size_t chan = diff.firstMismatch % 4u;
        const size_t x = pixel % static_cast<size_t>(exrW);
        const size_t y = pixel / static_cast<size_t>(exrW);
        std::printf("first mismatch: x=%zu y=%zu c=%zu exr=%.9g tile=%.9g\n",
                    x,
                    y,
                    chan,
                    exrRgba[diff.firstMismatch],
                    decoded->rgba[diff.firstMismatch]);
    }

    std::printf("verify result: %s\n", ok ? "PASS" : "FAIL");
    return ok ? 0 : 2;
}
