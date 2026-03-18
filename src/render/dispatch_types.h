#pragma once

#include "device/device_memory_view.h"

#include <cstdint>

namespace ybi
{

struct Scene;
struct ScenePool;
struct Camera;

enum class RenderKernelId : uint8_t
{
    PrimaryDiffuse = 0,
    AO = 1,
    PathTrace = 2,
};

struct DispatchParams
{
    ScenePool *scenePool = nullptr;
    Scene *rootScene = nullptr;
    const Camera *camera = nullptr;

    uint32_t width = 0;
    uint32_t height = 0;
    uint32_t spp = 1;
    uint64_t launchParamsDevice = 0;
    uint64_t launchParamsSize = 0;

    DeviceMemoryView<uint8_t> outputRGBA8 = {};
    DeviceMemoryView<float> outputLinearRgb = {};
    DeviceMemoryView<uint8_t> feedbackBytes = {};

    uint32_t frameIndex = 0;
    bool enableVirtualTexture = false;
};

} // namespace ybi
