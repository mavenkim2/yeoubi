#pragma once

#include "pxr/usd/ar/resolverContextBinder.h"
#include "pxr/usd/sdf/propertySpec.h"
#include "pxr/usd/usd/resolveInfo.h"
#include "pxr/usd/usdShade/shader.h"
#include "pxr/usd/usdShade/udimUtils.h"
#include "material/material.h"
#include <limits>

namespace ybi
{

inline pxr::SdfLayerHandle FindAssetPathAuthoringLayer(const pxr::UsdAttribute &sourceAttribute,
                                                       const pxr::SdfAssetPath &assetPath,
                                                       pxr::UsdTimeCode timeCode)
{
    const pxr::UsdResolveInfo resolveInfo = sourceAttribute.GetResolveInfo(timeCode);
    if (resolveInfo.GetSource() != pxr::UsdResolveInfoSourceDefault)
    {
        YBI_LOGFATAL("Asset paths with time-varying values are unsupported\n");
        return {};
    }

    for (const pxr::SdfPropertySpecHandle &propertySpec :
         sourceAttribute.GetPropertyStack(timeCode))
    {
        const pxr::SdfAttributeSpecHandle attributeSpec =
            pxr::TfDynamic_cast<pxr::SdfAttributeSpecHandle>(propertySpec);
        if (!attributeSpec || !attributeSpec->HasDefaultValue())
        {
            continue;
        }

        const pxr::VtValue &authoredValue = attributeSpec->GetDefaultValue();
        if (authoredValue.IsHolding<pxr::SdfAssetPath>() &&
            authoredValue.UncheckedGet<pxr::SdfAssetPath>().GetAssetPath() ==
                assetPath.GetAssetPath())
        {
            return attributeSpec->GetLayer();
        }
    }

    YBI_LOGFATAL("Could not find the layer that authored asset path: %s\n",
                 assetPath.GetAssetPath().c_str());
    return {};
}

inline void ParseUdimPaths(const pxr::SdfAssetPath &sdfAssetPath,
                           const pxr::UsdAttribute &sourceAttribute,
                           pxr::UsdTimeCode timeCode,
                           UdimTextureVector &filePaths)
{
    std::string rawPathStr = sdfAssetPath.GetAssetPath();
    std::string resolved = sdfAssetPath.GetResolvedPath();

    if (pxr::UsdShadeUdimUtils::IsUdimIdentifier(rawPathStr))
    {
        const pxr::SdfLayerHandle anchorLayer =
            FindAssetPathAuthoringLayer(sourceAttribute, sdfAssetPath, timeCode);
        if (!anchorLayer)
        {
            return;
        }

        const pxr::ArResolverContextBinder resolverContext(
            sourceAttribute.GetPrim().GetStage()->GetPathResolverContext());
        const auto tiles =
            pxr::UsdShadeUdimUtils::ResolveUdimTilePaths(rawPathStr, anchorLayer);

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
