#pragma once

#include "core/affine_transform.h"
#include "core/matrix4x4.h"
#include "device/device_texture.h"
#include "render/virtual_texture_interop.h"
#include "scene/device_scene.h"

namespace ybi
{

struct KernelParams
{
    // Camera
    Vec3f cameraPosition;
    float fovRadians;
    Matrix4x4 cameraFromRaster;
    AffineTransform renderFromCamera;
    float lensRadius;
    int imageWidth;
    int imageHeight;

    // BVH
    uint64_t bvhHandle;

    uint64_t imageHandle;

    // Geometry
    DeviceScene scene;
    int maxDepth;

    // Virtual textures
    DeviceTextureDescriptor *textureDescriptors;
    PageDescriptor *pageDescriptors;
    uint32_t *pageDescriptorBits;
    FeedbackRequest *feedbackRequests;
};

} // namespace ybi
