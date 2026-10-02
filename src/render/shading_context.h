#pragma once

#include "device/cuda_compat.h"
#include "device/device_texture.h"
#include "render/virtual_texture_interop.h"
#include "scene/device_scene.h"

namespace ybi
{

struct ShadingContext
{
    const DeviceScene &scene;

    const DeviceTextureDescriptor *textureDescriptors;
    const PageDescriptor *pageDescriptors;
    uint32_t *pageDescriptorBits;
    FeedbackRequest *feedbackRequests;

    uint64_t opaqueDielectricEnergyComplementTexture;
};

} // namespace ybi
