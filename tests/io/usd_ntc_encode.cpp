#include <pxr/usd/sdf/assetPath.h>
#include <pxr/usd/usd/primRange.h>
#include <pxr/usd/usd/stage.h>
#include <pxr/usd/usdShade/material.h>
#include <pxr/usd/usdShade/materialBindingAPI.h>
#include <pxr/usd/usdShade/shader.h>

#include <algorithm>
#include <cctype>
#include <cstdio>
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
    std::string ntcCli = "ntc-cli";
    float bitsPerPixel = 4.0f;
    bool noEncode = false;
};

struct ChannelTexture
{
    std::string texturePath;
    std::string swizzle;
};

struct MaterialChannels
{
    std::string materialPath;
    std::unordered_map<std::string, ChannelTexture> channels;
};

void PrintUsage(const char *exe)
{
    std::fprintf(stderr,
                 "Usage: %s <entry.usd[a|c]> <out_dir> [--ntc-cli <path>] [--bits-per-pixel <bpp>] [--no-encode]\n"
                 "Writes one NTC JSON manifest per material and invokes ntc-cli per material.\n",
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
        if (arg == "--bits-per-pixel")
        {
            if (i + 1 >= argc)
            {
                std::fprintf(stderr, "Missing value for --bits-per-pixel\n");
                return false;
            }
            out.bitsPerPixel = std::strtof(argv[++i], nullptr);
            if (!(out.bitsPerPixel > 0.0f))
            {
                std::fprintf(stderr, "Invalid --bits-per-pixel value\n");
                return false;
            }
            continue;
        }
        if (arg == "--no-encode")
        {
            out.noEncode = true;
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

std::string OutputNameToSwizzle(const pxr::TfToken &sourceName)
{
    const std::string name = sourceName.GetString();
    if (name == "r")
    {
        return "R";
    }
    if (name == "g")
    {
        return "G";
    }
    if (name == "b")
    {
        return "B";
    }
    if (name == "a")
    {
        return "A";
    }
    if (name == "rgb")
    {
        return "RGB";
    }
    if (name == "rgba")
    {
        return "RGBA";
    }
    return "";
}

std::string JsonEscape(const std::string &s)
{
    std::string out;
    out.reserve(s.size() + 8);
    for (char c : s)
    {
        if (c == '\\')
        {
            out += "\\\\";
        }
        else if (c == '"')
        {
            out += "\\\"";
        }
        else if (c == '\n')
        {
            out += "\\n";
        }
        else if (c == '\r')
        {
            out += "\\r";
        }
        else if (c == '\t')
        {
            out += "\\t";
        }
        else
        {
            out += c;
        }
    }
    return out;
}

std::string ToPortablePath(const std::string &path)
{
    return fs::path(path).generic_string();
}

std::string SemanticLabelForInput(const std::string &inputName)
{
    if (inputName == "diffuseColor")
    {
        return "Albedo";
    }
    if (inputName == "emissiveColor")
    {
        return "Emissive";
    }
    if (inputName == "normal")
    {
        return "Normal";
    }
    if (inputName == "roughness")
    {
        return "Roughness";
    }
    if (inputName == "metallic")
    {
        return "Metalness";
    }
    if (inputName == "occlusion")
    {
        return "Occlusion";
    }
    if (inputName == "opacity")
    {
        return "Alpha";
    }
    if (inputName == "specularColor")
    {
        return "SpecularColor";
    }
    return "";
}

std::string SemanticRangeForInput(const std::string &inputName, const std::string &swizzle)
{
    if (!swizzle.empty())
    {
        return swizzle;
    }
    if (inputName == "roughness" || inputName == "metallic" || inputName == "occlusion")
    {
        return "R";
    }
    if (inputName == "opacity")
    {
        return "A";
    }
    return "RGB";
}

bool IsSrgbInput(const std::string &inputName)
{
    return inputName == "diffuseColor" || inputName == "emissiveColor";
}

bool TryGetUVTextureFile(const pxr::UsdShadeInput &input,
                         ChannelTexture &out,
                         std::string &reasonOut)
{
    out = {};
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

    out.swizzle = OutputNameToSwizzle(sources[0].sourceName);

    const std::string path = AssetPathToString(fileAsset);
    if (path.empty())
    {
        reasonOut = "empty texture path";
        return false;
    }

    out.texturePath = path;
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

            ChannelTexture texture = {};
            std::string reason;
            if (TryGetUVTextureFile(input, texture, reason))
            {
                item.channels.emplace(input.GetBaseName().GetString(), std::move(texture));
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

    std::vector<std::pair<std::string, ChannelTexture>> sorted(mat.channels.begin(), mat.channels.end());
    std::sort(sorted.begin(), sorted.end(), [](const auto &a, const auto &b) {
        return a.first < b.first;
    });

    out << "{\n";
    out << "  \"textures\": [\n";
    for (size_t i = 0; i < sorted.size(); ++i)
    {
        const auto &kv = sorted[i];
        const std::string semantic = SemanticLabelForInput(kv.first);
        const std::string range = SemanticRangeForInput(kv.first, kv.second.swizzle);

        out << "    {\n";
        out << "      \"fileName\": \"" << JsonEscape(ToPortablePath(kv.second.texturePath)) << "\",\n";
        out << "      \"name\": \"" << JsonEscape(kv.first) << "\",\n";
        if (!kv.second.swizzle.empty())
        {
            out << "      \"channelSwizzle\": \"" << kv.second.swizzle << "\",\n";
        }
        out << "      \"isSRGB\": " << (IsSrgbInput(kv.first) ? "true" : "false");
        if (!semantic.empty())
        {
            out << ",\n";
            out << "      \"semantics\": {\n";
            out << "        \"" << semantic << "\": \"" << range << "\"\n";
            out << "      }\n";
        }
        else
        {
            out << "\n";
        }
        out << "    }";
        if (i + 1 < sorted.size())
        {
            out << ",";
        }
        out << "\n";
    }
    out << "  ]\n";
    out << "}\n";
    return true;
}

int RunNtc(const std::string &ntcCli, const fs::path &manifestPath, const fs::path &outPath, float bitsPerPixel)
{
    char bppBuf[64];
    std::snprintf(bppBuf, sizeof(bppBuf), "%.3f", bitsPerPixel);
    std::string cmd = "\"" + ntcCli + "\" --loadManifest \"" + manifestPath.string() +
                      "\" --generateMips --compress --bitsPerPixel " + bppBuf +
                      " --decompress --saveCompressed \"" + outPath.string() + "\"";
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

    if (!cli.noEncode)
    {
        std::printf("NTC encode executable: %s\n", cli.ntcCli.c_str());
    }

    for (const MaterialChannels &mat : materials)
    {
        if (mat.channels.empty())
        {
            continue;
        }
        ++withChannels;

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
            const int rc = RunNtc(cli.ntcCli, manifestPath, ntcOutPath, cli.bitsPerPixel);
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
            std::printf("  dry-run (--no-encode)\n");
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
