#pragma once

#include "core/matrix4x4.h"
#include "core/affine_transform.h"

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
    int maxDepth;
};

}
