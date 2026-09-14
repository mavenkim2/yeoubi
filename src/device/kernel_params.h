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
    bool pathTrace;
    uint32_t samplesPerPixel;
    uint32_t sampleOffset;
    Vec3f *accumulation;
    uint32_t accumulatedSamples;

    // OpenPBR lookup tables
    uint64_t opaqueDielectricEnergyComplementTexture;
    uint64_t opaqueDielectricAverageEnergyComplementTexture;

    // Virtual textures
    DeviceTextureDescriptor *textureDescriptors;
    PageDescriptor *pageDescriptors;
    uint32_t *pageDescriptorBits;
    FeedbackRequest *feedbackRequests;
};

} // namespace ybi
