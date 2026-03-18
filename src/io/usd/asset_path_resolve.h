#pragma once

#include <pxr/usd/sdf/assetPath.h>
#include <pxr/usd/sdf/layer.h>
#include <pxr/usd/sdf/layerUtils.h>
#include <pxr/usd/usd/attribute.h>
#include <pxr/usd/usdShade/connectableAPI.h>
#include <pxr/usd/usdShade/input.h>
#include <pxr/usd/usdShade/types.h>

#include <string>

namespace ybi
{

inline std::string ResolveUsdAssetPath(const pxr::UsdAttribute &attr,
                                       const pxr::SdfAssetPath &assetPath)
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

    const auto stack = attr.GetPropertyStack();
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

inline bool ResolveUsdShadeInputAssetPath(const pxr::UsdShadeInput &input,
                                          std::string *outPath,
                                          int depth = 0)
{
    constexpr int kMaxDepth = 16;

    if (!outPath)
    {
        return false;
    }
    outPath->clear();
    if (!input || depth > kMaxDepth)
    {
        return false;
    }

    pxr::SdfAssetPath assetPath;
    if (input.Get(&assetPath))
    {
        *outPath = ResolveUsdAssetPath(input.GetAttr(), assetPath);
        return !outPath->empty();
    }

    const auto sources = input.GetConnectedSources();
    if (sources.size() != 1)
    {
        return false;
    }

    const pxr::UsdShadeConnectionSourceInfo &source = sources[0];
    if (source.sourceType != pxr::UsdShadeAttributeType::Input)
    {
        return false;
    }

    const pxr::UsdShadeInput upstreamInput = source.source.GetInput(source.sourceName);
    if (!upstreamInput)
    {
        return false;
    }

    return ResolveUsdShadeInputAssetPath(upstreamInput, outPath, depth + 1);
}

} // namespace ybi
