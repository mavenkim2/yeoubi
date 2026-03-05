#include "tile_binary.h"
#include "texture/exr_io.h"
#include "texture/udim_utils.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <unordered_map>
#include <vector>

namespace
{

static constexpr uint32_t kUdimMin = 1001;
static constexpr uint32_t kUdimMax = 1100;

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
        if (!ybi::usd_ntc::TryFindUdimDigits(exrPath, udim, digitPos))
        {
            udim = kUdimMin;
        }
    }

    if (tilesBinPath.empty())
    {
        if (tilesDir.empty())
        {
            std::printf("missing tiles input: set --tiles-bin or --tiles-dir\n");
            return 1;
        }
        tilesBinPath = tilesDir + "/" + Sanitize(ybi::usd_ntc::StripUdimFromPath(exrPath)) + ".tiles.bin";
    }

    int exrW = 0;
    int exrH = 0;
    std::vector<float> exrRgba;
    std::string error;
    if (!ybi::texture::LoadExrRgba(exrPath, &exrW, &exrH, &exrRgba, &error, true))
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
