#pragma once

#include "pxr/usd/usdShade/shader.h"
#include "pxr/usd/usdShade/udimUtils.h"
#include "scene/material.h"

namespace ybi
{

inline void ParseUdimPaths(const pxr::SdfAssetPath &sdfAssetPath,
                           const pxr::UsdShadeShader &shader,
                           const pxr::UsdShadeInput &input,
                           UdimTextureVector &filePaths)
{
    std::string rawPathStr = sdfAssetPath.GetAssetPath();
    std::string resolved = sdfAssetPath.GetResolvedPath();

    if (pxr::UsdShadeUdimUtils::IsUdimIdentifier(rawPathStr))
    {
        pxr::SdfLayerHandle anchorLayer;
        pxr::SdfPropertySpecHandleVector stack = input.GetAttr().GetPropertyStack();

        if (!stack.empty())
        {
            anchorLayer = stack.front()->GetLayer();
        }
        else
        {
            anchorLayer = shader.GetPrim().GetStage()->GetRootLayer();
        }

        auto tiles = pxr::UsdShadeUdimUtils::ResolveUdimTilePaths(rawPathStr, anchorLayer);

        for (const auto &[resolvedPath, tile] : tiles)
        {
            auto value = std::stoul(tile);
            YBI_ASSERT(value <= std::numeric_limits<uint32_t>::max());

            filePaths.emplace_back(resolvedPath, static_cast<uint32_t>(value));
        }
    }
    else if (!resolved.empty())
    {
        filePaths.emplace_back(resolved, ~0u);
    }
    else
    {
        YBI_LOGFATAL("Failure when getting paths from pxr::SdfAssetPath\n");
    }
}

} // namespace ybi
