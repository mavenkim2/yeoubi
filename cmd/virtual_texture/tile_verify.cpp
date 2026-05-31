#include "exr_mips.h"
#include "tile_binary.h"
#include "tile_pixels.h"
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
    std::printf("  verifies the full mip chain, not just mip 0\n");
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
        if (!ybi::texture::TryFindUdimDigits(exrPath, udim, digitPos))
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
        tilesBinPath = tilesDir + "/" + Sanitize(ybi::texture::StripUdimFromPath(exrPath)) + ".tiles.bin";
    }

    std::string error;
    ybi::tileprep::ExrMipChainResult sourceMipChain = {};
    if (!ybi::tileprep::LoadExrMipChain(exrPath, true, &sourceMipChain, &error))
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
    if (sourceMipChain.mipLevels.empty())
    {
        std::printf("verify failed: source mip chain is empty for %s\n", exrPath.c_str());
        return 2;
    }
    if (decoded->mipLevels.size() != sourceMipChain.mipLevels.size())
    {
        std::printf("verify failed: mip count mismatch exr=%zu tile=%zu udim=%u\n",
                    sourceMipChain.mipLevels.size(),
                    decoded->mipLevels.size(),
                    udim);
        return 2;
    }

    std::printf("verify stats: exr=%s tiles=%s udim=%u\n", exrPath.c_str(), tilesBinPath.c_str(), udim);
    std::printf("verify stats: width=%u height=%u mipCount=%zu channels=%u storedMips=%s\n",
                sourceMipChain.mipLevels[0].width,
                sourceMipChain.mipLevels[0].height,
                sourceMipChain.mipLevels.size(),
                ybi::TextureFormatChannelCount(decoded->pixelFormat),
                sourceMipChain.hasStoredMipLevels ? "yes" : "no");

    bool allOk = true;
    for (size_t mipIndex = 0; mipIndex < sourceMipChain.mipLevels.size(); ++mipIndex)
    {
        const ybi::tilebin::UdimMipImage &srcMip = sourceMipChain.mipLevels[mipIndex];
        const ybi::tilebin::UdimMipImage &dstMip = decoded->mipLevels[mipIndex];
        if (srcMip.level != dstMip.level || srcMip.width != dstMip.width ||
            srcMip.height != dstMip.height || dstMip.pixelFormat != decoded->pixelFormat)
        {
            std::printf("verify mip %zu: dimension mismatch exr=%ux%u tile=%ux%u levels=%u/%u\n",
                        mipIndex,
                        srcMip.width,
                        srcMip.height,
                        dstMip.width,
                        dstMip.height,
                        srcMip.level,
                        dstMip.level);
            allOk = false;
            continue;
        }

        std::vector<float> srcPixels;
        if (!ybi::tileprep::ConvertPixelsToPixelFormat(
                decoded->pixelFormat, srcMip.rgba, &srcPixels, &error))
        {
            std::printf("verify failed: %s\n", error.c_str());
            return 2;
        }

        ybi::tilebin::DiffStats diff = {};
        const bool mipOk = ybi::tilebin::DiffImagesExact(srcPixels, dstMip.rgba, eps, &diff);
        std::printf("verify mip %zu: width=%u height=%u maxAbs=%.9g meanAbs=%.9g rmse=%.9g mismatches=%llu/%zu %s\n",
                    mipIndex,
                    srcMip.width,
                    srcMip.height,
                    diff.maxAbs,
                    diff.meanAbs,
                    diff.rmse,
                    static_cast<unsigned long long>(diff.mismatchCount),
                    srcPixels.size(),
                    mipOk ? "PASS" : "FAIL");
        if (!mipOk && diff.firstMismatch < srcPixels.size())
        {
            const size_t channels =
                static_cast<size_t>(ybi::TextureFormatChannelCount(
                    decoded->pixelFormat));
            const size_t pixel = diff.firstMismatch / channels;
            const size_t chan = diff.firstMismatch % channels;
            const size_t x = pixel % static_cast<size_t>(srcMip.width);
            const size_t y = pixel / static_cast<size_t>(srcMip.width);
            std::printf("first mismatch mip %zu: x=%zu y=%zu c=%zu exr=%.9g tile=%.9g\n",
                        mipIndex,
                        x,
                        y,
                        chan,
                        srcPixels[diff.firstMismatch],
                        dstMip.rgba[diff.firstMismatch]);
        }
        allOk = allOk && mipOk;
    }

    std::printf("verify result: %s\n", allOk ? "PASS" : "FAIL");
    return allOk ? 0 : 2;
}
