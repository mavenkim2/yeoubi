#pragma once

#include "device/device_memory.h"
#include <memory>

namespace ybi
{

struct Device;
struct DeviceKernel;
struct WavefrontParams;

class WavefrontIntegrator
{
public:
    explicit WavefrontIntegrator(Device *device);
    ~WavefrontIntegrator();

    bool RenderSamples(const WavefrontParams &params);

private:
    Device *device;
    std::unique_ptr<DeviceKernel> countKernels;
    std::unique_ptr<DeviceKernel> selectSampleIndices;
    DeviceBuffer kernelCountsBuffer;
    DeviceBuffer numSampleIndicesBuffer;
};

} // namespace ybi
