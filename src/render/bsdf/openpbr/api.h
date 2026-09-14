#pragma once

#include "render/bsdf/openpbr/fuzz.h"

namespace ybi
{

YBI_GPU bool SampleOpenPBR(const OpenPbrBsdf &bsdf,
                           const Vec3f &worldViewDirection,
                           const Vec3f &throughput,
                           Vec2f u,
                           Vec3f &worldLightDirection,
                           Vec3f &weight)
{
    const Vec3f localViewDirection = bsdf.tangentFrame.WorldToLocal(worldViewDirection);
    Vec3f localLightDirection;
    const bool success =
        SampleFuzzLayer(bsdf, localViewDirection, throughput, u, localLightDirection, weight);
    if (success)
    {
        worldLightDirection = bsdf.tangentFrame.LocalToWorld(localLightDirection);
    }
    return success;
}

YBI_GPU Vec3f EvaluateOpenPBR(const OpenPbrBsdf &bsdf,
                              const Vec3f &worldViewDirection,
                              const Vec3f &worldLightDirection)
{
    const Vec3f localViewDirection = bsdf.tangentFrame.WorldToLocal(worldViewDirection);
    const Vec3f localLightDirection = bsdf.tangentFrame.WorldToLocal(worldLightDirection);
    return EvaluateOpenPBRLocal(bsdf, localViewDirection, localLightDirection);
}

YBI_GPU float PDFOpenPBR(const OpenPbrBsdf &bsdf,
                         const Vec3f &worldViewDirection,
                         const Vec3f &worldLightDirection)
{
    const Vec3f localViewDirection = bsdf.tangentFrame.WorldToLocal(worldViewDirection);
    const Vec3f localLightDirection = bsdf.tangentFrame.WorldToLocal(worldLightDirection);
    return PDFOpenPBRLocal(bsdf, localViewDirection, localLightDirection);
}

} // namespace ybi
