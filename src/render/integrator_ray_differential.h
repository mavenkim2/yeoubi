#pragma once

#include "render/integrator_common.h"
#include "render/launch_params.h"

namespace ybi
{
namespace render
{
namespace integrator
{

struct RayDifferential
{
    Vec3 origin = Vec3(0.0f, 0.0f, 0.0f);
    Vec3 dir = Vec3(0.0f, 0.0f, 1.0f);
    Vec3 originX = Vec3(0.0f, 0.0f, 0.0f);
    Vec3 originY = Vec3(0.0f, 0.0f, 0.0f);
    Vec3 dirX = Vec3(0.0f, 0.0f, 1.0f);
    Vec3 dirY = Vec3(0.0f, 0.0f, 1.0f);
    bool valid = false;
};

YBI_INTEGRATOR_HD Vec3 ComputePerspectiveCameraDirection(const LaunchParams &params,
                                                         float rasterX,
                                                         float rasterY,
                                                         unsigned int width,
                                                         unsigned int height)
{
    if (width == 0u || height == 0u)
    {
        return Vec3(0.0f, 0.0f, 1.0f);
    }

    const float invWidth = 1.0f / static_cast<float>(width);
    const float invHeight = 1.0f / static_cast<float>(height);
    const float ndcX = rasterX * invWidth * 2.0f - 1.0f;
    const float ndcY = 1.0f - rasterY * invHeight * 2.0f;
    return Normalize(Vec3(params.cameraU.x * ndcX + params.cameraV.x * ndcY + params.cameraW.x,
                          params.cameraU.y * ndcX + params.cameraV.y * ndcY + params.cameraW.y,
                          params.cameraU.z * ndcX + params.cameraV.z * ndcY + params.cameraW.z));
}

YBI_INTEGRATOR_HD RayDifferential InitPerspectiveRayDifferential(const LaunchParams &params,
                                                                 float pixelX,
                                                                 float pixelY,
                                                                 unsigned int width,
                                                                 unsigned int height)
{
    RayDifferential rayDiff = {};
    const Vec3 cameraOrigin(params.cameraOrigin.x, params.cameraOrigin.y, params.cameraOrigin.z);
    rayDiff.origin = cameraOrigin;
    rayDiff.originX = cameraOrigin;
    rayDiff.originY = cameraOrigin;

    if (width == 0u || height == 0u)
    {
        return rayDiff;
    }

    rayDiff.dir =
        ComputePerspectiveCameraDirection(params, pixelX + 0.5f, pixelY + 0.5f, width, height);
    rayDiff.dirX =
        ComputePerspectiveCameraDirection(params, pixelX + 1.5f, pixelY + 0.5f, width, height);
    rayDiff.dirY =
        ComputePerspectiveCameraDirection(params, pixelX + 0.5f, pixelY + 1.5f, width, height);
    rayDiff.valid = true;
    return rayDiff;
}

} // namespace integrator
} // namespace render
} // namespace ybi
