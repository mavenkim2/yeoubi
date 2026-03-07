#pragma once

#include "render/integrator_common.h"
#include "render/launch_params.h"

namespace ybi
{
namespace render
{
namespace integrator
{

YBI_INTEGRATOR_HD int ComputeTextureMipCount(int width, int height)
{
    int maxDim = MaxInt(width, height);
    maxDim = MaxInt(maxDim, 1);
    int mipCount = 1;
    while (maxDim > 1)
    {
        maxDim = MaxInt(maxDim >> 1, 1);
        mipCount++;
    }
    return mipCount;
}

YBI_INTEGRATOR_HD bool ProjectWorldToRaster(const LaunchParams &params,
                                            const Vec3 &worldPoint,
                                            Vec2 *outRaster)
{
    if (!outRaster || params.width <= 0 || params.height <= 0)
    {
        return false;
    }

    const Vec3 cameraOrigin =
        MakeVec3(params.cameraOrigin.x, params.cameraOrigin.y, params.cameraOrigin.z);
    const Vec3 cameraU = MakeVec3(params.cameraU.x, params.cameraU.y, params.cameraU.z);
    const Vec3 cameraV = MakeVec3(params.cameraV.x, params.cameraV.y, params.cameraV.z);
    const Vec3 cameraW = MakeVec3(params.cameraW.x, params.cameraW.y, params.cameraW.z);
    const Vec3 rel = Sub(worldPoint, cameraOrigin);

    const float det = Dot(cameraU, Cross(cameraV, cameraW));
    if (fabsf(det) <= 1e-8f)
    {
        return false;
    }

    const float alpha = Dot(rel, Cross(cameraV, cameraW)) / det;
    const float beta = Dot(cameraU, Cross(rel, cameraW)) / det;
    const float gamma = Dot(cameraU, Cross(cameraV, rel)) / det;
    if (fabsf(gamma) <= 1e-8f)
    {
        return false;
    }

    const float ndcX = alpha / gamma;
    const float ndcY = beta / gamma;
    outRaster->x = (ndcX + 1.0f) * 0.5f * float(params.width) - 0.5f;
    outRaster->y = (1.0f - ndcY) * 0.5f * float(params.height) - 0.5f;
    return true;
}

YBI_INTEGRATOR_HD bool TryComputeTextureMipLevel(const LaunchParams &params,
                                                 const HitInfo &hit,
                                                 const UV2 &uv0,
                                                 const UV2 &uv1,
                                                 const UV2 &uv2,
                                                 int textureWidth,
                                                 int textureHeight,
                                                 unsigned int *outMip)
{
    if (!outMip || textureWidth <= 0 || textureHeight <= 0 || !hit.hasWorldTriangle)
    {
        return false;
    }

    Vec2 p0 = {};
    Vec2 p1 = {};
    Vec2 p2 = {};
    if (!ProjectWorldToRaster(params, hit.worldTri0, &p0) ||
        !ProjectWorldToRaster(params, hit.worldTri1, &p1) ||
        !ProjectWorldToRaster(params, hit.worldTri2, &p2))
    {
        return false;
    }

    const float dx1 = p1.x - p0.x;
    const float dy1 = p1.y - p0.y;
    const float dx2 = p2.x - p0.x;
    const float dy2 = p2.y - p0.y;
    const float det = dx1 * dy2 - dy1 * dx2;
    if (fabsf(det) <= 1e-8f)
    {
        return false;
    }

    const float invDet = 1.0f / det;
    const float du1 = uv1.x - uv0.x;
    const float dv1 = uv1.y - uv0.y;
    const float du2 = uv2.x - uv0.x;
    const float dv2 = uv2.y - uv0.y;
    const float dudx = (du1 * dy2 - du2 * dy1) * invDet;
    const float dudy = (du2 * dx1 - du1 * dx2) * invDet;
    const float dvdx = (dv1 * dy2 - dv2 * dy1) * invDet;
    const float dvdy = (dv2 * dx1 - dv1 * dx2) * invDet;

    const float rhoX = sqrtf((dudx * float(textureWidth)) * (dudx * float(textureWidth)) +
                             (dvdx * float(textureHeight)) * (dvdx * float(textureHeight)));
    const float rhoY = sqrtf((dudy * float(textureWidth)) * (dudy * float(textureWidth)) +
                             (dvdy * float(textureHeight)) * (dvdy * float(textureHeight)));
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

    const int maxMip = MaxInt(ComputeTextureMipCount(textureWidth, textureHeight) - 1, 0);
    if (mip > maxMip)
    {
        mip = maxMip;
    }
    mip = ClampInt(mip, 0, 15);
    *outMip = static_cast<unsigned int>(mip);
    return true;
}

} // namespace integrator
} // namespace render
} // namespace ybi
