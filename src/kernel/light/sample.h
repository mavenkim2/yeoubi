#pragma once

#include "device/cuda_compat.h"
#include "kernel/light/envmap.h"
#include "kernel/light/rect.h"
#include "kernel/scene.h"

namespace ybi
{

YBI_GPU uint32_t GetLightCount(const DeviceScene &scene)
{
    return scene.numLights + scene.numEnvmapLights;
}

YBI_GPU bool SampleLight(const DeviceScene &scene,
                         const Vec3f &hitPosition,
                         Vec3f u,
                         LightSample &sample)
{
    const uint32_t numLights = GetLightCount(scene);
    YBI_ASSERT(numLights > 0);

    const uint32_t lightIndex =
        Min(static_cast<uint32_t>(u.x * numLights), numLights - 1u);

    if (lightIndex < scene.numLights)
    {
        const DeviceLight &light = scene.deviceLights[lightIndex];
        Vec3f lightPosition;
        sample.direction =
            SampleRectLight(light, hitPosition, Vec2f(u.y, u.z), lightPosition, sample.pdf);
        sample.radiance = light.radiance;
        sample.distance = Length(lightPosition - hitPosition);
    }
    else
    {
        const DeviceEnvmapLight &light = scene.envmapLights[lightIndex - scene.numLights];
        if (!SampleEnvmapLight(
                light, Vec2f(u.y, u.z), sample.direction, sample.radiance, sample.pdf))
            return false;
        sample.distance = FLOAT_INFINITY;
    }

    sample.pdf /= static_cast<float>(numLights);
    return sample.pdf != 0.f;
}

YBI_GPU float PDFLight(const DeviceScene &scene,
                       uint32_t lightIndex,
                       const Vec3f &hitPosition,
                       const Vec3f &lightPosition)
{
    YBI_ASSERT(lightIndex < scene.numLights);

    const DeviceLight &light = scene.deviceLights[lightIndex];
    return PDFRectLight(light, hitPosition, lightPosition) /
           static_cast<float>(GetLightCount(scene));
}

YBI_GPU float PDFEnvmapLight(const DeviceScene &scene,
                             uint32_t lightIndex,
                             const Vec3f &direction)
{
    return PDFEnvmapLight(scene.envmapLights[lightIndex], direction) /
           static_cast<float>(GetLightCount(scene));
}

} // namespace ybi
