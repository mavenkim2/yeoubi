#pragma once

#include "render/integrator_common.h"
#include "util/math_common.h"

namespace ybi
{
namespace render
{
namespace integrator
{

YBI_DEVICE int ComputeTextureMipCount(int width, int height)
{
    int maxDim = Max(width, height);
    maxDim = Max(maxDim, 1);
    int mipCount = 1;
    while (maxDim > 1)
    {
        maxDim = Max(maxDim >> 1, 1);
        mipCount++;
    }
    return mipCount;
}

YBI_DEVICE bool TryComputeTextureMipLevel(const HitInfo &hit,
                                          int textureWidth,
                                          int textureHeight,
                                          unsigned int *outMip)
{
    if (!outMip || textureWidth <= 0 || textureHeight <= 0 || !hit.hasTextureDifferentials)
    {
        return false;
    }

    const float rhoX =
        sqrtf((hit.dSdx * float(textureWidth)) * (hit.dSdx * float(textureWidth)) +
              (hit.dTdx * float(textureHeight)) * (hit.dTdx * float(textureHeight)));
    const float rhoY =
        sqrtf((hit.dSdy * float(textureWidth)) * (hit.dSdy * float(textureWidth)) +
              (hit.dTdy * float(textureHeight)) * (hit.dTdy * float(textureHeight)));
    const float rho = rhoX > rhoY ? rhoX : rhoY;
    if (!(rho >= 0.0f) || rho > 1.0e30f)
    {
        return false;
    }

    int mip = 0;
    if (rho > 1.0f)
    {
        const float lambda = log2f(rho);
        if (!(lambda >= 0.0f) || lambda > 1.0e30f)
        {
            return false;
        }
        mip = static_cast<int>(floorf(lambda));
    }

    const int maxMip = Max(ComputeTextureMipCount(textureWidth, textureHeight) - 1, 0);
    if (mip > maxMip)
    {
        mip = maxMip;
    }
    mip = Clamp(mip, 0, 15);
    *outMip = static_cast<unsigned int>(mip);
    return true;
}

} // namespace integrator
} // namespace render
} // namespace ybi
