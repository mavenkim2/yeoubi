#pragma once

#include "util/math_common.h"
#include <cstdint>

namespace ybi
{
namespace render
{

YBI_DEVICE float EvalAOVisibility(int visibleSamples, int totalSamples)
{
    const int denom = totalSamples > 0 ? totalSamples : 1;
    const int visible = visibleSamples < 0 ? 0 : (visibleSamples > denom ? denom : visibleSamples);
    return static_cast<float>(visible) / static_cast<float>(denom);
}

YBI_DEVICE void EvalPrimaryDiffuseCore(bool sampled,
                                       float sampledR,
                                       float sampledG,
                                       float sampledB,
                                       float &outR,
                                       float &outG,
                                       float &outB)
{
    if (!sampled)
    {
        outR = 0.0f;
        outG = 0.0f;
        outB = 0.0f;
        return;
    }
    outR = sampledR;
    outG = sampledG;
    outB = sampledB;
}

YBI_DEVICE uint32_t PackRGB8(float r, float g, float b)
{
    const uint32_t ru = static_cast<uint32_t>(Clamp(r, 0.0f, 1.0f) * 255.0f + 0.5f);
    const uint32_t gu = static_cast<uint32_t>(Clamp(g, 0.0f, 1.0f) * 255.0f + 0.5f);
    const uint32_t bu = static_cast<uint32_t>(Clamp(b, 0.0f, 1.0f) * 255.0f + 0.5f);
    return ru | (gu << 8) | (bu << 16);
}

} // namespace render

} // namespace ybi
