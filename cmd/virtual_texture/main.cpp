#include "shared.h"

#include <pxr/usd/usd/stage.h>

#include <chrono>
#include <cstdio>

int main(int argc, char **argv)
{
    const auto runStart = std::chrono::steady_clock::now();
    const auto LogTotalRuntime = [&runStart]() {
        const auto runEnd = std::chrono::steady_clock::now();
        const double ms = std::chrono::duration<double, std::milli>(runEnd - runStart).count();
        std::fprintf(stderr, "Tile stage: total runtime %.3f ms\n", ms);
    };

    Cli cli = {};
    if (!ParseCli(argc, argv, cli))
    {
        LogTotalRuntime();
        return 2;
    }

    std::fprintf(stderr, "Tile stage: open USD: %s\n", cli.usdPath.c_str());
    pxr::UsdStageRefPtr stage = pxr::UsdStage::Open(cli.usdPath);
    if (!stage)
    {
        std::fprintf(stderr, "Failed to open USD stage: %s\n", cli.usdPath.c_str());
        LogTotalRuntime();
        return 1;
    }
    std::fprintf(stderr, "Tile stage: collect material textures\n");
    const std::vector<MaterialChannels> materials = CollectMaterialChannels(stage, cli.purposes);
    std::fprintf(stderr, "Tile stage: collect material textures done (%zu materials)\n", materials.size());

    std::string tileError;
    if (!PrepareTexturesForStreamingTiles(materials, cli, &tileError))
    {
        std::fprintf(stderr, "Tile prep failed: %s\n", tileError.c_str());
        LogTotalRuntime();
        return 1;
    }

    LogTotalRuntime();
    return 0;
}
