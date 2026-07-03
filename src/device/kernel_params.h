#pragma once

#include "core/affine_transform.h"
#include "core/matrix4x4.h"
#include "scene/device_scene.h"

namespace ybi
{

struct KernelParams
{
    // Camera
    Matrix4x4 cameraFromRaster;
    AffineTransform renderFromCamera;
    float lensRadius;

    // BVH
    uint64_t bvhHandle;

    uint64_t imageHandle;

    // Geometry
    DeviceScene scene;
    int maxDepth;
};

} // namespace ybi
