#pragma once

#include <cstdint>

namespace ybi
{
namespace render
{

#if defined(__CUDACC__)
#define YBI_RENDER_HD __host__ __device__ inline
#else
#define YBI_RENDER_HD inline
#endif

YBI_RENDER_HD float Clamp01(float x)
{
    if (x < 0.0f)
    {
        return 0.0f;
    }
    if (x > 1.0f)
    {
        return 1.0f;
    }
    return x;
}

YBI_RENDER_HD float EvalAOVisibility(int visibleSamples, int totalSamples)
{
    const int denom = totalSamples > 0 ? totalSamples : 1;
    const int visible = visibleSamples < 0 ? 0 : (visibleSamples > denom ? denom : visibleSamples);
    return static_cast<float>(visible) / static_cast<float>(denom);
}

YBI_RENDER_HD void EvalPrimaryDiffuseCore(bool sampled,
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

YBI_RENDER_HD uint32_t PackRGB8(float r, float g, float b)
{
    const uint32_t ru = static_cast<uint32_t>(Clamp01(r) * 255.0f + 0.5f);
    const uint32_t gu = static_cast<uint32_t>(Clamp01(g) * 255.0f + 0.5f);
    const uint32_t bu = static_cast<uint32_t>(Clamp01(b) * 255.0f + 0.5f);
    return ru | (gu << 8) | (bu << 16);
}

#undef YBI_RENDER_HD

} // namespace render

} // namespace ybi
