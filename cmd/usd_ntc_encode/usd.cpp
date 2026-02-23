#include "shared.h"

#include <pxr/usd/sdf/assetPath.h>
#include <pxr/usd/sdf/layer.h>
#include <pxr/usd/sdf/layerUtils.h>
#include <pxr/usd/sdf/valueTypeName.h>
#include <pxr/usd/usd/primRange.h>
#include <pxr/usd/usdGeom/imageable.h>
#include <pxr/usd/usdGeom/mesh.h>
#include <pxr/usd/usdGeom/tokens.h>
#include <pxr/usd/usdShade/material.h>
#include <pxr/usd/usdShade/materialBindingAPI.h>
#include <pxr/usd/usdShade/shader.h>

#include <algorithm>
#include <cctype>
#include <cstdio>
#include <cstdlib>

namespace
{

bool IsConnected(const pxr::UsdShadeInput &input)
{
    return !input.GetConnectedSources().empty();
}

bool ParsePurposeList(const std::string &csv, std::vector<std::string> &outPurposes)
{
    outPurposes.clear();
    size_t begin = 0;
    while (begin <= csv.size())
    {
        const size_t comma = csv.find(',', begin);
        const size_t end = (comma == std::string::npos) ? csv.size() : comma;
        size_t b = begin;
        while (b < end && std::isspace((unsigned char)csv[b]))
        {
            ++b;
        }
        size_t e = end;
        while (e > b && std::isspace((unsigned char)csv[e - 1]))
        {
            --e;
        }
        std::string token = csv.substr(b, e - b);
        for (char &c : token)
        {
            c = (char)std::tolower((unsigned char)c);
        }
        if (!token.empty())
        {
            if (token != "default" && token != "render" && token != "proxy" &&
                token != "guide")
            {
                return false;
            }
            if (std::find(outPurposes.begin(), outPurposes.end(), token) == outPurposes.end())
            {
                outPurposes.push_back(token);
            }
        }
        if (comma == std::string::npos)
        {
            break;
        }
        begin = comma + 1;
    }
    return !outPurposes.empty();
}

bool IsPurposeAllowed(const pxr::TfToken &purpose, const std::vector<std::string> &allowedPurposes)
{
    const char *name = "";
    if (purpose == pxr::UsdGeomTokens->default_)
    {
        name = "default";
    }
    else if (purpose == pxr::UsdGeomTokens->render)
    {
        name = "render";
    }
    else if (purpose == pxr::UsdGeomTokens->proxy)
    {
        name = "proxy";
    }
    else if (purpose == pxr::UsdGeomTokens->guide)
    {
        name = "guide";
    }
    for (const std::string &allowed : allowedPurposes)
    {
        if (allowed == name)
        {
            return true;
        }
    }
    return false;
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
                 "Usage: %s <entry.usd[a|c]> <out_dir> [--out-usd path] [--bits-per-pixel <bpp>] [--training-steps <n>] [--steps-per-iter <n>] [--cuda-device <n>] [--max-materials <n>] [--purposes <csv>] [--no-encode]\n",
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
        if (arg == "--out-usd")
        {
            if (i + 1 >= argc)
            {
                std::fprintf(stderr, "Missing value for --out-usd\n");
                return false;
            }
            out.outUsdPath = argv[++i];
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
        if (arg == "--purposes")
        {
            if (i + 1 >= argc)
            {
                std::fprintf(stderr, "Missing value for --purposes\n");
                return false;
            }
            if (!ParsePurposeList(argv[++i], out.purposes))
            {
                std::fprintf(stderr,
                             "Invalid --purposes. expected comma-separated values from: default,render,proxy,guide\n");
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

std::vector<MaterialChannels> CollectMaterialChannels(const pxr::UsdStageRefPtr &stage,
                                                      const std::vector<std::string> &purposes)
{
    std::unordered_map<std::string, pxr::UsdShadeMaterial> uniqueMaterials;

    const pxr::Usd_PrimFlagsPredicate pred =
        pxr::UsdPrimIsActive && pxr::UsdPrimIsLoaded && !pxr::UsdPrimIsAbstract;
    for (const pxr::UsdPrim &prim : stage->GetPseudoRoot().GetFilteredDescendants(pred))
    {
        if (!prim.IsA<pxr::UsdGeomMesh>() || prim.IsInPrototype())
        {
            continue;
        }
        pxr::UsdGeomImageable imageable(prim);
        if (!imageable)
        {
            continue;
        }
        const pxr::TfToken purpose = imageable.ComputePurpose();
        if (!IsPurposeAllowed(purpose, purposes))
        {
            continue;
        }

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

bool WriteNtcBindingsToUsd(const pxr::UsdStageRefPtr &stage,
                           const std::unordered_map<std::string, std::string> &materialToNtcFile,
                           const std::string &outUsdPath,
                           std::string *outError)
{
    if (!stage)
    {
        if (outError)
        {
            *outError = "null USD stage";
        }
        return false;
    }
    if (outUsdPath.empty())
    {
        if (outError)
        {
            *outError = "empty output USD path";
        }
        return false;
    }

    const pxr::TfToken fileToken("ybi:ntc:diffuseFile");
    const pxr::TfToken texNameToken("ybi:ntc:diffuseTextureName");
    int authoredCount = 0;
    for (const auto &kv : materialToNtcFile)
    {
        const std::string &materialPath = kv.first;
        const std::string &ntcFile = kv.second;
        if (materialPath.empty() || ntcFile.empty())
        {
            continue;
        }

        pxr::UsdPrim materialPrim = stage->GetPrimAtPath(pxr::SdfPath(materialPath));
        if (!materialPrim || !materialPrim.IsA<pxr::UsdShadeMaterial>())
        {
            std::printf("NTC USD write: material missing %s\n", materialPath.c_str());
            continue;
        }

        pxr::UsdAttribute fileAttr = materialPrim.CreateAttribute(
            fileToken, pxr::SdfValueTypeNames->Asset, true);
        fileAttr.Set(pxr::SdfAssetPath(ntcFile));

        pxr::UsdAttribute texNameAttr = materialPrim.CreateAttribute(
            texNameToken, pxr::SdfValueTypeNames->String, true);
        texNameAttr.Set(std::string("diffuseColor"));
        authoredCount++;
    }

    if (!stage->Export(outUsdPath))
    {
        if (outError)
        {
            *outError = "stage export failed: " + outUsdPath;
        }
        return false;
    }

std::printf("NTC USD write: authored bindings=%d exported=%s\n",
                authoredCount,
                outUsdPath.c_str());
    return true;
}
