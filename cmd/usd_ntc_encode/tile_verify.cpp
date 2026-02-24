#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <string>
#include <type_traits>
#include <vector>

#define TINYEXR_IMPLEMENTATION
#include <tinyexr.h>

namespace
{

struct TileFileHeader
{
    char magic[8];
    uint32_t version = 1;
    uint32_t tileSize = 0;
    uint32_t imageWidth = 0;
    uint32_t imageHeight = 0;
    uint32_t channels = 4;
    uint32_t tileCountX = 0;
    uint32_t tileCountY = 0;
    uint32_t tileCount = 0;
    uint32_t elementType = 1; // 1 = float32
};

struct TileRecord
{
    uint32_t tileX = 0;
    uint32_t tileY = 0;
    uint32_t width = 0;
    uint32_t height = 0;
    uint64_t byteOffset = 0;
    uint64_t byteSize = 0;
};

static_assert(std::is_trivially_copyable<TileFileHeader>::value, "TileFileHeader must be POD");
static_assert(std::is_trivially_copyable<TileRecord>::value, "TileRecord must be POD");

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
        out = "material";
    }
    return out;
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

bool ReadTileImage(const std::string &tilePath, int &width, int &height, std::vector<float> &rgba, std::string &error)
{
    std::ifstream in(tilePath, std::ios::binary);
    if (!in.is_open())
    {
        error = "failed to open tile file: " + tilePath;
        return false;
    }

    TileFileHeader header = {};
    in.read(reinterpret_cast<char *>(&header), sizeof(header));
    if (!in.good())
    {
        error = "failed reading tile header: " + tilePath;
        return false;
    }
    if (std::memcmp(header.magic, "YBITILE1", 8) != 0)
    {
        error = "bad tile magic: " + tilePath;
        return false;
    }
    if (header.channels != 4 || header.elementType != 1 || header.tileSize == 0)
    {
        error = "unsupported tile format in: " + tilePath;
        return false;
    }
    if (header.tileCount != header.tileCountX * header.tileCountY)
    {
        error = "tile count mismatch in: " + tilePath;
        return false;
    }

    std::vector<TileRecord> records(header.tileCount);
    in.read(reinterpret_cast<char *>(records.data()), static_cast<std::streamsize>(records.size() * sizeof(TileRecord)));
    if (!in.good())
    {
        error = "failed reading tile records: " + tilePath;
        return false;
    }

    width = static_cast<int>(header.imageWidth);
    height = static_cast<int>(header.imageHeight);
    rgba.assign(static_cast<size_t>(width) * static_cast<size_t>(height) * 4u, 0.0f);

    for (const TileRecord &r : records)
    {
        if (r.tileX >= header.tileCountX || r.tileY >= header.tileCountY)
        {
            error = "tile record out of range in: " + tilePath;
            return false;
        }
        if (r.width == 0 || r.height == 0 || r.width > header.tileSize || r.height > header.tileSize)
        {
            error = "tile record size invalid in: " + tilePath;
            return false;
        }

        const uint64_t expectedBytes = static_cast<uint64_t>(r.width) * static_cast<uint64_t>(r.height) * 4u * sizeof(float);
        if (r.byteSize != expectedBytes)
        {
            error = "tile byte size mismatch in: " + tilePath;
            return false;
        }

        std::vector<float> tile(static_cast<size_t>(r.width) * static_cast<size_t>(r.height) * 4u);
        in.seekg(static_cast<std::streamoff>(r.byteOffset), std::ios::beg);
        in.read(reinterpret_cast<char *>(tile.data()), static_cast<std::streamsize>(r.byteSize));
        if (!in.good())
        {
            error = "failed reading tile payload in: " + tilePath;
            return false;
        }

        const uint32_t x0 = r.tileX * header.tileSize;
        const uint32_t y0 = r.tileY * header.tileSize;
        if (x0 + r.width > header.imageWidth || y0 + r.height > header.imageHeight)
        {
            error = "tile bounds overflow in: " + tilePath;
            return false;
        }

        for (uint32_t y = 0; y < r.height; ++y)
        {
            const float *src = tile.data() + static_cast<size_t>(y) * static_cast<size_t>(r.width) * 4u;
            float *dst = rgba.data() +
                         (static_cast<size_t>(y0 + y) * static_cast<size_t>(header.imageWidth) + static_cast<size_t>(x0)) * 4u;
            std::memcpy(dst, src, static_cast<size_t>(r.width) * 4u * sizeof(float));
        }
    }

    return true;
}

void PrintUsage(const char *exe)
{
    std::printf("Usage:\n");
    std::printf("  %s --exr <path.exr> --tiles-bin <path.tiles.bin> [--eps <float>]\n", exe);
    std::printf("  %s --exr <path.exr> --tiles-dir <dir> [--eps <float>]\n", exe);
    std::printf("Notes:\n");
    std::printf("  --tiles-dir auto-resolves to <dir>/<sanitize(exr)>.tiles.bin\n");
}

} // namespace

int main(int argc, char **argv)
{
    std::string exrPath;
    std::string tilesBinPath;
    std::string tilesDir;
    float eps = 1e-7f;

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

    if (tilesBinPath.empty())
    {
        if (tilesDir.empty())
        {
            std::printf("missing tiles input: set --tiles-bin or --tiles-dir\n");
            return 1;
        }
        tilesBinPath = tilesDir + "/" + Sanitize(exrPath) + ".tiles.bin";
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

    int tileW = 0;
    int tileH = 0;
    std::vector<float> tileRgba;
    if (!ReadTileImage(tilesBinPath, tileW, tileH, tileRgba, error))
    {
        std::printf("verify failed: %s\n", error.c_str());
        return 2;
    }

    if (exrW != tileW || exrH != tileH || exrRgba.size() != tileRgba.size())
    {
        std::printf("verify failed: dimension mismatch exr=%dx%d tile=%dx%d\n", exrW, exrH, tileW, tileH);
        return 2;
    }

    double sumSq = 0.0;
    double sumAbs = 0.0;
    float maxAbs = 0.0f;
    uint64_t mismatchCount = 0;
    size_t firstMismatch = exrRgba.size();
    for (size_t i = 0; i < exrRgba.size(); ++i)
    {
        const float diff = std::fabs(exrRgba[i] - tileRgba[i]);
        maxAbs = std::max(maxAbs, diff);
        sumAbs += diff;
        sumSq += static_cast<double>(diff) * static_cast<double>(diff);
        if (diff > eps)
        {
            ++mismatchCount;
            if (firstMismatch == exrRgba.size())
            {
                firstMismatch = i;
            }
        }
    }

    const double denom = exrRgba.empty() ? 1.0 : static_cast<double>(exrRgba.size());
    const double meanAbs = sumAbs / denom;
    const double rmse = std::sqrt(sumSq / denom);
    std::printf("verify stats: exr=%s tiles=%s\n", exrPath.c_str(), tilesBinPath.c_str());
    std::printf("verify stats: width=%d height=%d channels=4\n", exrW, exrH);
    std::printf("verify stats: eps=%g maxAbs=%.9g meanAbs=%.9g rmse=%.9g mismatches=%llu/%zu\n",
                eps,
                maxAbs,
                meanAbs,
                rmse,
                static_cast<unsigned long long>(mismatchCount),
                exrRgba.size());

    if (firstMismatch < exrRgba.size())
    {
        const size_t pixel = firstMismatch / 4u;
        const size_t chan = firstMismatch % 4u;
        const size_t x = pixel % static_cast<size_t>(exrW);
        const size_t y = pixel / static_cast<size_t>(exrW);
        std::printf("first mismatch: x=%zu y=%zu c=%zu exr=%.9g tile=%.9g\n",
                    x,
                    y,
                    chan,
                    exrRgba[firstMismatch],
                    tileRgba[firstMismatch]);
    }

    if (mismatchCount > 0)
    {
        std::printf("verify result: FAIL\n");
        return 2;
    }

    std::printf("verify result: PASS\n");
    return 0;
}
