#include "shared.h"

#include <pxr/usd/usd/stage.h>

#include <libntc/ntc.h>
#include <libntc/wrappers.h>

#include <cstdio>
#include <cstdint>
#include <filesystem>
#include <string>
#include <unordered_map>
#include <vector>

namespace
{

void PrintProgress(int done, int total)
{
    const int safeTotal = total > 0 ? total : 1;
    const int width = 30;
    const int filled = (done * width) / safeTotal;
    const int remaining = total - done;

    std::printf("[");
    for (int i = 0; i < width; ++i)
    {
        std::printf("%c", i < filled ? '=' : ' ');
    }
    std::printf("] %d/%d texture sets, %d remaining\n", done, total, remaining < 0 ? 0 : remaining);
}

double ToMiB(uint64_t bytes)
{
    return static_cast<double>(bytes) / (1024.0 * 1024.0);
}

} // namespace

int main(int argc, char **argv)
{
    Cli cli = {};
    if (!ParseCli(argc, argv, cli))
    {
        return 2;
    }

    std::fprintf(stderr, "NTC stage: open USD: %s\n", cli.usdPath.c_str());
    pxr::UsdStageRefPtr stage = pxr::UsdStage::Open(cli.usdPath);
    if (!stage)
    {
        std::fprintf(stderr, "Failed to open USD stage: %s\n", cli.usdPath.c_str());
        return 1;
    }
    std::fprintf(stderr, "NTC stage: open USD done\n");

    std::filesystem::create_directories(cli.outDir);
    if (cli.outUsdPath.empty())
    {
        const fs::path inPath(cli.usdPath);
        const std::string stem = inPath.stem().string();
        cli.outUsdPath = (fs::path(cli.outDir) / (stem + ".ntc.usda")).string();
    }
    fs::path outUsdDir = fs::path(cli.outUsdPath).parent_path();
    if (!outUsdDir.empty())
    {
        std::filesystem::create_directories(outUsdDir);
    }

    std::fprintf(stderr, "NTC stage: collect material textures\n");
    std::fprintf(stderr, "NTC stage: purposes=");
    for (size_t i = 0; i < cli.purposes.size(); ++i)
    {
        std::fprintf(stderr, "%s%s", i ? "," : "", cli.purposes[i].c_str());
    }
    std::fprintf(stderr, "\n");
    const std::vector<MaterialChannels> materials = CollectMaterialChannels(stage, cli.purposes);
    std::fprintf(stderr, "NTC stage: collect material textures done (%zu materials)\n", materials.size());

    if (cli.prepareTiles)
    {
        std::string tileError;
        if (!PrepareTexturesForStreamingTiles(materials, cli, &tileError))
        {
            std::fprintf(stderr, "Tile prep failed: %s\n", tileError.c_str());
            return 1;
        }
        return 0;
    }

    ntc::ContextWrapper context;
    if (!cli.noEncode)
    {
        std::fprintf(stderr, "NTC stage: create NTC context (cudaDevice=%d)\n", cli.cudaDevice);
        ntc::ContextParameters contextParams = {};
        contextParams.cudaDevice = cli.cudaDevice;
        const ntc::Status status = ntc::CreateContext(context.ptr(), contextParams);
        if (status != ntc::Status::Ok)
        {
            std::fprintf(stderr,
                         "Failed to create NTC context: %s (%s)\n",
                         ntc::StatusToString(status),
                         ntc::GetLastErrorMessage());
            return 1;
        }
        std::fprintf(stderr, "NTC stage: create NTC context done\n");
    }

    int encodedCount = 0;
    int encodeFailCount = 0;
    std::unordered_map<std::string, std::string> materialToNtcFile;
    int withChannels = 0;
    uint64_t totalSourceFileBytes = 0;
    uint64_t totalDecodedBytes = 0;
    uint64_t totalNtcBytes = 0;
    int totalToProcess = 0;
    for (const MaterialChannels &mat : materials)
    {
        if (mat.channels.empty())
        {
            continue;
        }
        if (cli.maxMaterials > 0 && totalToProcess >= cli.maxMaterials)
        {
            break;
        }
        ++totalToProcess;
    }
    int processed = 0;

    PrintProgress(0, totalToProcess);

    for (const MaterialChannels &mat : materials)
    {
        if (mat.channels.empty())
        {
            continue;
        }
        if (cli.maxMaterials > 0 && withChannels >= cli.maxMaterials)
        {
            break;
        }

        ++withChannels;
        ++processed;

        const std::string base = Sanitize(mat.materialPath);
        const fs::path ntcOutPath = fs::path(cli.outDir) / (base + ".ntc");
        std::printf("material %s channels=%zu\n", mat.materialPath.c_str(), mat.channels.size());

        if (!cli.noEncode)
        {
            float actualBpp = 0.0f;
            std::string reason;
            EncodeStats encodeStats = {};
            if (EncodeMaterial(context,
                               mat,
                               ntcOutPath,
                               cli.bitsPerPixel,
                               cli.trainingSteps,
                               cli.stepsPerIter,
                               processed,
                               totalToProcess,
                               encodeStats,
                               actualBpp,
                               reason))
            {
                ++encodedCount;
                materialToNtcFile[mat.materialPath] = ntcOutPath.filename().string();
                totalSourceFileBytes += encodeStats.sourceFileBytes;
                totalDecodedBytes += encodeStats.decodedBytes;
                totalNtcBytes += encodeStats.ntcBytes;
                std::printf("  ntc encoded -> %s (bpp=%.3f)\n", ntcOutPath.string().c_str(), actualBpp);
                std::printf("  footprint source=%.2f MiB decoded=%.2f MiB ntc=%.2f MiB\n",
                            ToMiB(encodeStats.sourceFileBytes),
                            ToMiB(encodeStats.decodedBytes),
                            ToMiB(encodeStats.ntcBytes));
            }
            else
            {
                ++encodeFailCount;
                totalSourceFileBytes += encodeStats.sourceFileBytes;
                totalDecodedBytes += encodeStats.decodedBytes;
                std::printf("  ntc encode failed: %s\n", reason.c_str());
                std::printf("  footprint source=%.2f MiB decoded=%.2f MiB ntc=0.00 MiB\n",
                            ToMiB(encodeStats.sourceFileBytes),
                            ToMiB(encodeStats.decodedBytes));
            }
        }
        else
        {
            std::printf("  dry-run (--no-encode)\n");
        }

        PrintProgress(processed, totalToProcess);
    }

    std::printf("NTC summary materials=%zu materialsWithChannels=%d encoded=%d encodeFailed=%d\n",
                materials.size(),
                withChannels,
                encodedCount,
                encodeFailCount);
    if (!cli.noEncode)
    {
        std::printf("NTC footprint summary source=%.2f MiB decoded=%.2f MiB ntc=%.2f MiB\n",
                    ToMiB(totalSourceFileBytes),
                    ToMiB(totalDecodedBytes),
                    ToMiB(totalNtcBytes));
        if (!materialToNtcFile.empty())
        {
            std::string writeError;
            if (!WriteNtcBindingsToUsd(stage, materialToNtcFile, cli.outUsdPath, &writeError))
            {
                std::fprintf(stderr, "Failed to write NTC bindings to USD: %s\n", writeError.c_str());
                return 1;
            }
        }
        else
        {
            std::printf("NTC USD write: skipped (no encoded material bindings)\n");
        }
    }
    else
    {
        std::printf("NTC USD write: skipped in --no-encode mode\n");
    }

    return 0;
}
