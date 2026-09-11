#pragma once

#include "device/cuda_compat.h"
#include "device/kernel_params.h"
#include "render/lights/rect.h"

namespace ybi
{

YBI_GPU Vec3f SampleLight(const KernelParams &params,
                          const Vec3f &hitPosition,
                          Vec3f u,
                          int &sampledLightIndex,
                          Vec3f &lightPosition,
                          float &pdf)
{
    YBI_ASSERT(params.scene.numLights > 0);

    const uint32_t lightIndex =
        Min(static_cast<uint32_t>(u.x * params.scene.numLights), params.scene.numLights - 1u);

    const DeviceLight &light = params.scene.deviceLights[lightIndex];

    sampledLightIndex = static_cast<int>(lightIndex);
    const Vec3f lightDirection =
        SampleRectLight(light, hitPosition, Vec2f(u.y, u.z), lightPosition, pdf);
    pdf /= static_cast<float>(params.scene.numLights);
    return lightDirection;
}

YBI_GPU float PDFLight(const KernelParams &params,
                       uint32_t lightIndex,
                       const Vec3f &hitPosition,
                       const Vec3f &lightPosition)
{
    YBI_ASSERT(lightIndex < params.scene.numLights);

    const DeviceLight &light = params.scene.deviceLights[lightIndex];
    return PDFRectLight(light, hitPosition, lightPosition) /
           static_cast<float>(params.scene.numLights);
}

} // namespace ybi
