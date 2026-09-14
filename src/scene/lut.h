#pragma once

#include "device/device_texture.h"

namespace ybi
{

struct Device;

DeviceTexture
LoadOpenPbrOpaqueDielectricEnergyComplementTexture(Device &device, const char *executablePath);

DeviceTexture LoadOpenPbrOpaqueDielectricAverageEnergyComplementTexture(
    Device &device, const char *executablePath);

} // namespace ybi
