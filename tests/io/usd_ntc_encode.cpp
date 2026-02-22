#include <pxr/usd/sdf/assetPath.h>
#include <pxr/usd/usd/primRange.h>
#include <pxr/usd/usd/stage.h>
#include <pxr/usd/usdShade/material.h>
#include <pxr/usd/usdShade/materialBindingAPI.h>
#include <pxr/usd/usdShade/shader.h>

#include <algorithm>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <string>
#include <unordered_map>
#include <vector>

namespace fs = std::filesystem;

namespace
{

struct Cli
{
    std::string usdPath;
    std::string outDir;
    std::string ntcCli;
};

struct MaterialChannels
{
    std::string materialPath;
    std::unordered_map<std::string, std::string> channels;
};

void PrintUsage(const char *exe)
{
    std::fprintf(stderr,
                 "Usage: %s <entry.usd[a|c]> <out_dir> [--ntc-cli <path>]\n"
                 "Writes one manifest per material and optionally runs NTC CLI per material.\n",
                 exe);
}

bool ParseCli(int argc, char **argv, Cli &out)
{
    if (argc < 3)
    {
        PrintUsage(argv[0]);
        return false;
    }
    out.usdPath = argv[1];
    out.outDir = argv[2];

    for (int i = 3; i < argc; ++i)
    {
        const std::string arg = argv[i];
        if (arg == "--ntc-cli")
        {
            if (i + 1 >= argc)
            {
                std::fprintf(stderr, "Missing value for --ntc-cli\n");
                return false;
            }
            out.ntcCli = argv[++i];
            continue;
        }
        std::fprintf(stderr, "Unknown option: %s\n", arg.c_str());
        return false;
    }
    return true;
}

std::string AssetPathToString(const pxr::SdfAssetPath &asset)
{
    std::string resolved = asset.GetResolvedPath();
    if (!resolved.empty())
    {
        return resolved;
    }
    return asset.GetAssetPath();
}

bool IsConnected(const pxr::UsdShadeInput &input)
{
    return !input.GetConnectedSources().empty();
}

bool TryGetUVTextureFile(const pxr::UsdShadeInput &input,
                         std::string &pathOut,
                         std::string &reasonOut)
{
    pathOut.clear();
    reasonOut.clear();

    auto sources = input.GetConnectedSources();
    if (sources.size() != 1)
    {
        reasonOut = "unsupported connection count";
        return false;
    }

    pxr::UsdPrim sourcePrim = sources[0].source.GetPrim();
    if (!sourcePrim.IsA<pxr::UsdShadeShader>())
    {
        reasonOut = "source is not shader";
        return false;
    }

    pxr::UsdShadeShader textureShader(sourcePrim);
    pxr::TfToken shaderId;
    textureShader.GetShaderId(&shaderId);
    if (shaderId != pxr::TfToken("UsdUVTexture"))
    {
        // TODO: handle non-UsdUVTexture image nodes (e.g. MaterialX image nodes).
        reasonOut = "unsupported texture shader: " + shaderId.GetString();
        return false;
    }

    pxr::UsdShadeInput fileInput = textureShader.GetInput(pxr::TfToken("file"));
    if (!fileInput)
    {
        // TODO: handle texture nodes with alternate file-bearing input names.
        reasonOut = "UsdUVTexture missing input:file";
        return false;
    }

    pxr::SdfAssetPath fileAsset;
    if (!fileInput.Get(&fileAsset))
    {
        // TODO: handle file input connections or non-asset value types.
        reasonOut = "input:file has no readable SdfAssetPath";
        return false;
    }

    const std::string path = AssetPathToString(fileAsset);
    if (path.empty())
    {
        reasonOut = "empty texture path";
        return false;
    }

    pathOut = path;
    return true;
}

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

std::vector<MaterialChannels> CollectMaterialChannels(const pxr::UsdStageRefPtr &stage)
{
    std::unordered_map<std::string, pxr::UsdShadeMaterial> uniqueMaterials;

    const pxr::Usd_PrimFlagsPredicate pred =
        pxr::UsdPrimIsActive && pxr::UsdPrimIsLoaded && !pxr::UsdPrimIsAbstract;
    for (const pxr::UsdPrim &prim : stage->GetPseudoRoot().GetFilteredDescendants(pred))
    {
        pxr::UsdShadeMaterialBindingAPI bindingApi(prim);
        pxr::UsdShadeMaterial material = bindingApi.ComputeBoundMaterial(pxr::UsdShadeTokens->full);
        if (!material)
        {
            continue;
        }
        const std::string materialPath = material.GetPath().GetString();
        uniqueMaterials.emplace(materialPath, material);
    }

    std::vector<MaterialChannels> out;
    out.reserve(uniqueMaterials.size());

    for (const auto &it : uniqueMaterials)
    {
        MaterialChannels item = {};
        item.materialPath = it.first;

        pxr::UsdShadeMaterial material = it.second;
        pxr::UsdShadeShader surface = material.ComputeSurfaceSource();
        if (!surface)
        {
            std::printf("NTC gather: material %s has no surface source\n", item.materialPath.c_str());
            out.push_back(std::move(item));
            continue;
        }

        pxr::TfToken surfaceId;
        surface.GetShaderId(&surfaceId);
        if (surfaceId != pxr::TfToken("UsdPreviewSurface"))
        {
            // TODO: support non-UsdPreviewSurface material graphs.
            std::printf("NTC gather: material %s unsupported surface shader %s\n",
                        item.materialPath.c_str(),
                        surfaceId.GetText());
            out.push_back(std::move(item));
            continue;
        }

        for (const pxr::UsdShadeInput &input : surface.GetInputs())
        {
            if (!IsConnected(input))
            {
                continue;
            }

            std::string texturePath;
            std::string reason;
            if (TryGetUVTextureFile(input, texturePath, reason))
            {
                item.channels.emplace(input.GetBaseName().GetString(), texturePath);
            }
            else
            {
                std::printf("NTC gather: material %s input %s skipped (%s)\n",
                            item.materialPath.c_str(),
                            input.GetBaseName().GetText(),
                            reason.c_str());
            }
        }

        out.push_back(std::move(item));
    }

    std::sort(out.begin(), out.end(), [](const MaterialChannels &a, const MaterialChannels &b) {
        return a.materialPath < b.materialPath;
    });
    return out;
}

bool WriteManifest(const fs::path &path, const MaterialChannels &mat)
{
    std::ofstream out(path);
    if (!out.is_open())
    {
        return false;
    }

    out << "material " << mat.materialPath << "\n";
    out << "channel_count " << mat.channels.size() << "\n";

    std::vector<std::pair<std::string, std::string>> sorted(mat.channels.begin(), mat.channels.end());
    std::sort(sorted.begin(), sorted.end());
    for (const auto &kv : sorted)
    {
        out << "channel " << kv.first << " " << kv.second << "\n";
    }
    return true;
}

int RunNtc(const std::string &ntcCli, const fs::path &manifestPath, const fs::path &outPath)
{
    std::string cmd = "\"" + ntcCli + "\" encode --manifest \"" + manifestPath.string() +
                      "\" --output \"" + outPath.string() + "\"";
    return std::system(cmd.c_str());
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

    fs::create_directories(cli.outDir);

    const std::vector<MaterialChannels> materials = CollectMaterialChannels(stage);

    int manifestCount = 0;
    int encodedCount = 0;
    int encodeFailCount = 0;
    int withChannels = 0;

    for (const MaterialChannels &mat : materials)
    {
        if (mat.channels.empty())
        {
            continue;
        }
        ++withChannels;

        const std::string base = Sanitize(mat.materialPath);
        const fs::path manifestPath = fs::path(cli.outDir) / (base + ".ntc_manifest.txt");
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

        if (!cli.ntcCli.empty())
        {
            const int rc = RunNtc(cli.ntcCli, manifestPath, ntcOutPath);
            if (rc == 0)
            {
                ++encodedCount;
                std::printf("  ntc encoded -> %s\n", ntcOutPath.string().c_str());
            }
            else
            {
                ++encodeFailCount;
                std::printf("  ntc encode failed rc=%d for %s\n", rc, mat.materialPath.c_str());
            }
        }
        else
        {
            std::printf("  dry-run (set --ntc-cli to encode)\n");
        }
    }

    std::printf("NTC summary materials=%zu materialsWithChannels=%d manifests=%d encoded=%d encodeFailed=%d\n",
                materials.size(),
                withChannels,
                manifestCount,
                encodedCount,
                encodeFailCount);

    return 0;
}
