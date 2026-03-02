#include "../usd_ntc_encode/shared.h"

#include <pxr/usd/usd/stage.h>

#include <cstdio>

int main(int argc, char **argv)
{
    Cli cli = {};
    if (!ParseCli(argc, argv, cli))
    {
        return 2;
    }

    cli.prepareTiles = true;
    cli.noEncode = true;

    std::fprintf(stderr, "Tile stage: open USD: %s\n", cli.usdPath.c_str());
    pxr::UsdStageRefPtr stage = pxr::UsdStage::Open(cli.usdPath);
    if (!stage)
    {
        std::fprintf(stderr, "Failed to open USD stage: %s\n", cli.usdPath.c_str());
        return 1;
    }
    std::fprintf(stderr, "Tile stage: collect material textures\n");
    const std::vector<MaterialChannels> materials = CollectMaterialChannels(stage, cli.purposes);
    std::fprintf(stderr, "Tile stage: collect material textures done (%zu materials)\n", materials.size());

    std::string tileError;
    if (!PrepareTexturesForStreamingTiles(materials, cli, &tileError))
    {
        std::fprintf(stderr, "Tile prep failed: %s\n", tileError.c_str());
        return 1;
    }

    return 0;
}
