#pragma once

#include "device/device_texture.h"

namespace ybi
{

struct Device;

DeviceTexture
LoadOpenPbrOpaqueDielectricEnergyComplementTexture(Device &device, const char *executablePath);

} // namespace ybi
