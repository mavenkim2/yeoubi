#include "shared.h"

#include <pxr/usd/sdf/assetPath.h>
#include <pxr/usd/sdf/layer.h>
#include <pxr/usd/sdf/layerUtils.h>
#include <pxr/usd/usd/primRange.h>
#include <pxr/usd/usdShade/material.h>
#include <pxr/usd/usdShade/materialBindingAPI.h>
#include <pxr/usd/usdShade/shader.h>

#include <algorithm>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <fstream>

namespace
{

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

std::string ResolveAssetPath(const pxr::UsdShadeInput &fileInput, const pxr::SdfAssetPath &assetPath)
{
    if (!assetPath.GetResolvedPath().empty())
    {
        return assetPath.GetResolvedPath();
    }

    const std::string raw = assetPath.GetAssetPath();
    if (raw.empty())
    {
        return {};
    }

    const auto stack = fileInput.GetAttr().GetPropertyStack();
    for (const auto &spec : stack)
    {
        if (!spec)
        {
            continue;
        }
        pxr::SdfLayerHandle layer = spec->GetLayer();
        if (!layer)
        {
            continue;
        }

        std::string resolved = pxr::SdfResolveAssetPathRelativeToLayer(layer, raw);
        if (!resolved.empty())
        {
            return resolved;
        }

        std::string anchored = pxr::SdfComputeAssetPathRelativeToLayer(layer, raw);
        if (!anchored.empty())
        {
            return anchored;
        }
    }

    return raw;
}

bool IsSrgbInput(const std::string &inputName)
{
    return inputName == "diffuseColor" || inputName == "emissiveColor";
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
        reasonOut = "unsupported texture shader: " + shaderId.GetString();
        return false;
    }

    pxr::UsdShadeInput fileInput = textureShader.GetInput(pxr::TfToken("file"));
    if (!fileInput)
    {
        reasonOut = "UsdUVTexture missing input:file";
        return false;
    }

    pxr::SdfAssetPath fileAsset;
    if (!fileInput.Get(&fileAsset))
    {
        reasonOut = "input:file has no readable SdfAssetPath";
        return false;
    }

    const std::string path = ResolveAssetPath(fileInput, fileAsset);
    if (path.empty())
    {
        reasonOut = "empty texture path";
        return false;
    }

    out.texturePath = path;
    out.swizzle = OutputNameToSwizzle(sources[0].sourceName);
    return true;
}

} // namespace

void PrintUsage(const char *exe)
{
    std::fprintf(stderr,
                 "Usage: %s <entry.usd[a|c]> <out_dir> [--bits-per-pixel <bpp>] [--training-steps <n>] [--steps-per-iter <n>] [--cuda-device <n>] [--max-materials <n>] [--no-encode]\n",
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
                std::fprintf(stderr, "Invalid --bits-per-pixel\n");
                return false;
            }
            continue;
        }
        if (arg == "--training-steps")
        {
            if (i + 1 >= argc)
            {
                std::fprintf(stderr, "Missing value for --training-steps\n");
                return false;
            }
            out.trainingSteps = std::atoi(argv[++i]);
            if (out.trainingSteps <= 0)
            {
                std::fprintf(stderr, "Invalid --training-steps\n");
                return false;
            }
            continue;
        }
        if (arg == "--steps-per-iter")
        {
            if (i + 1 >= argc)
            {
                std::fprintf(stderr, "Missing value for --steps-per-iter\n");
                return false;
            }
            out.stepsPerIter = std::atoi(argv[++i]);
            if (out.stepsPerIter <= 0)
            {
                std::fprintf(stderr, "Invalid --steps-per-iter\n");
                return false;
            }
            continue;
        }
        if (arg == "--cuda-device")
        {
            if (i + 1 >= argc)
            {
                std::fprintf(stderr, "Missing value for --cuda-device\n");
                return false;
            }
            out.cudaDevice = std::atoi(argv[++i]);
            continue;
        }
        if (arg == "--max-materials")
        {
            if (i + 1 >= argc)
            {
                std::fprintf(stderr, "Missing value for --max-materials\n");
                return false;
            }
            out.maxMaterials = std::atoi(argv[++i]);
            if (out.maxMaterials < 0)
            {
                std::fprintf(stderr, "Invalid --max-materials\n");
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

            ChannelTexture texture;
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
