#pragma once

#include "core/affine_transform.h"
#include "core/matrix4x4.h"
#include "device/device_texture.h"
#include "render/types.h"
#include "render/virtual_texture_interop.h"
#include "scene/device_scene.h"

namespace ybi
{

struct WavefrontParams
{
    // Camera
    Vec3f cameraPosition;
    float fovRadians;
    Matrix4x4 cameraFromRaster;
    AffineTransform renderFromCamera;
    float lensRadius;
    int imageWidth;
    int imageHeight;
    uint32_t pixelOffsetX;
    uint32_t pixelOffsetY;

    // BVH
    uint64_t bvhHandle;

    uint64_t imageHandle;
    Vec3f *kernelOutput;

    // Geometry
    DeviceScene scene;
    int maxDepth;
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

    // Wavefront
    SampleState *sampleStates;
    uint32_t sampleIndex;
    uint32_t maxNumSampleState;

    IntersectionSampleState *intersectionInfos;
    uint32_t numIntersectionInfos;
};

} // namespace ybi
