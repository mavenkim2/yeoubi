#include "shared.h"

#include <pxr/usd/usd/stage.h>

#include <libntc/ntc.h>
#include <libntc/wrappers.h>

#include <cstdio>
#include <filesystem>
#include <string>
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

} // namespace

int main(int argc, char **argv)
{
    Cli cli = {};
    if (!ParseCli(argc, argv, cli))
    {
        return 2;
    }

    pxr::UsdStageRefPtr stage = pxr::UsdStage::Open(cli.usdPath);
    if (!stage)
    {
        std::fprintf(stderr, "Failed to open USD stage: %s\n", cli.usdPath.c_str());
        return 1;
    }

    std::filesystem::create_directories(cli.outDir);

    const std::vector<MaterialChannels> materials = CollectMaterialChannels(stage);

    ntc::ContextWrapper context;
    if (!cli.noEncode)
    {
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
    }

    int manifestCount = 0;
    int encodedCount = 0;
    int encodeFailCount = 0;
    int withChannels = 0;
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
        const fs::path manifestPath = fs::path(cli.outDir) / (base + ".ntc_manifest.json");
        const fs::path ntcOutPath = fs::path(cli.outDir) / (base + ".ntc");

        if (!WriteManifest(manifestPath, mat))
        {
            std::fprintf(stderr, "Failed to write manifest: %s\n", manifestPath.string().c_str());
            return 1;
        }
        ++manifestCount;

        std::printf("material %s channels=%zu manifest=%s\n",
                    mat.materialPath.c_str(),
                    mat.channels.size(),
                    manifestPath.string().c_str());

        if (!cli.noEncode)
        {
            float actualBpp = 0.0f;
            std::string reason;
            if (EncodeMaterial(context,
                               mat,
                               ntcOutPath,
                               cli.bitsPerPixel,
                               cli.trainingSteps,
                               cli.stepsPerIter,
                               processed,
                               totalToProcess,
                               actualBpp,
                               reason))
            {
                ++encodedCount;
                std::printf("  ntc encoded -> %s (bpp=%.3f)\n", ntcOutPath.string().c_str(), actualBpp);
            }
            else
            {
                ++encodeFailCount;
                std::printf("  ntc encode failed: %s\n", reason.c_str());
            }
        }
        else
        {
            std::printf("  dry-run (--no-encode)\n");
        }

        PrintProgress(processed, totalToProcess);
    }

    std::printf("NTC summary materials=%zu materialsWithChannels=%d manifests=%d encoded=%d encodeFailed=%d\n",
                materials.size(),
                withChannels,
                manifestCount,
                encodedCount,
                encodeFailCount);

    return 0;
}
