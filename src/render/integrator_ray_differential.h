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
                                                         float rasterY)
{
    const Vec3 cameraPoint =
        TransformPointPerspective(params.cameraFromRaster, Vec3(rasterX, rasterY, 0.0f));
    if (!(cameraPoint.x == cameraPoint.x && cameraPoint.y == cameraPoint.y &&
          cameraPoint.z == cameraPoint.z))
    {
        return Vec3(0.0f, 0.0f, 1.0f);
    }

    const Vec3 worldPoint = TransformPointAffine(params.worldFromCamera, cameraPoint);
    const Vec3 worldOrigin =
        Vec3(params.cameraOrigin.x, params.cameraOrigin.y, params.cameraOrigin.z);
    return Normalize(worldPoint - worldOrigin);
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

    rayDiff.dir = ComputePerspectiveCameraDirection(params, pixelX + 0.5f, pixelY + 0.5f);
    rayDiff.dirX = ComputePerspectiveCameraDirection(params, pixelX + 1.5f, pixelY + 0.5f);
    rayDiff.dirY = ComputePerspectiveCameraDirection(params, pixelX + 0.5f, pixelY + 1.5f);
    rayDiff.valid = true;
    return rayDiff;
}

} // namespace integrator
} // namespace render
} // namespace ybi
