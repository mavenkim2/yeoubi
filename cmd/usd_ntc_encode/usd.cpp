#include "shared.h"

#include <pxr/usd/sdf/assetPath.h>
#include <pxr/usd/sdf/layer.h>
#include <pxr/usd/sdf/layerUtils.h>
#include <pxr/usd/sdf/valueTypeName.h>
#include <pxr/usd/usd/primRange.h>
#include <pxr/usd/usdGeom/gprim.h>
#include <pxr/usd/usdGeom/imageable.h>
#include <pxr/usd/usdGeom/mesh.h>
#include <pxr/usd/usdGeom/tokens.h>
#include <pxr/usd/usdShade/material.h>
#include <pxr/usd/usdShade/materialBindingAPI.h>
#include <pxr/usd/usdShade/shader.h>

#include <algorithm>
#include <cctype>
#include <cstdio>
#include <unordered_set>

namespace
{

bool IsConnected(const pxr::UsdShadeInput &input)
{
    return !input.GetConnectedSources().empty();
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

bool IsPurposeAllowedForPrim(const pxr::UsdPrim &prim, const std::vector<std::string> &allowedPurposes)
{
    pxr::UsdGeomImageable imageable(prim);
    if (!imageable)
    {
        // Non-imageable instance roots can still bind to prototypes that contain
        // renderable gprims; keep them to avoid dropping reachable materials.
        return true;
    }
    return IsPurposeAllowed(imageable.ComputePurpose(), allowedPurposes);
}

void AddBoundMaterialForPrim(const pxr::UsdPrim &prim,
                             std::unordered_map<std::string, pxr::UsdShadeMaterial> &materials)
{
    pxr::UsdShadeMaterialBindingAPI bindingApi(prim);
    pxr::UsdShadeMaterial material = bindingApi.ComputeBoundMaterial(pxr::UsdShadeTokens->full);
    if (!material)
    {
        return;
    }
    const std::string materialPath = material.GetPath().GetString();
    materials.emplace(materialPath, material);
}

} // namespace

std::vector<MaterialChannels> CollectMaterialChannels(const pxr::UsdStageRefPtr &stage,
                                                      const std::vector<std::string> &purposes)
{
    std::unordered_map<std::string, pxr::UsdShadeMaterial> uniqueMaterials;
    std::vector<pxr::UsdPrim> prototypeQueue;
    std::unordered_set<std::string> queuedPrototypePaths;
    prototypeQueue.reserve(64);
    queuedPrototypePaths.reserve(64);

    const pxr::Usd_PrimFlagsPredicate pred =
        pxr::UsdPrimIsActive && pxr::UsdPrimIsLoaded && !pxr::UsdPrimIsAbstract;

    const auto QueuePrototype = [&](const pxr::UsdPrim &instancePrim) {
        if (!instancePrim || !instancePrim.IsInstance())
        {
            return;
        }
        const pxr::UsdPrim prototype = instancePrim.GetPrototype();
        if (!prototype)
        {
            return;
        }
        const std::string prototypePath = prototype.GetPath().GetString();
        if (queuedPrototypePaths.emplace(prototypePath).second)
        {
            prototypeQueue.push_back(prototype);
        }
    };

    for (const pxr::UsdPrim &prim : stage->GetPseudoRoot().GetFilteredDescendants(pred))
    {
        if (prim.IsInPrototype())
        {
            continue;
        }

        if (prim.IsA<pxr::UsdGeomGprim>() && IsPurposeAllowedForPrim(prim, purposes))
        {
            AddBoundMaterialForPrim(prim, uniqueMaterials);
        }

        if (prim.IsInstance() && IsPurposeAllowedForPrim(prim, purposes))
        {
            QueuePrototype(prim);
        }
    }

    for (size_t i = 0; i < prototypeQueue.size(); ++i)
    {
        for (const pxr::UsdPrim &prim : pxr::UsdPrimRange(prototypeQueue[i], pred))
        {
            if (prim.IsA<pxr::UsdGeomGprim>())
            {
                AddBoundMaterialForPrim(prim, uniqueMaterials);
            }
            if (prim.IsInstance())
            {
                QueuePrototype(prim);
            }
        }
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
