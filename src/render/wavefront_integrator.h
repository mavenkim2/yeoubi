#pragma once

#include "device/device_memory.h"
#include "kernel/state/types.h"
#include <memory>
#include <vector>

namespace ybi
{

struct Device;
struct DeviceKernel;
struct KernelParams;

class WavefrontIntegrator
{
public:
    explicit WavefrontIntegrator(Device *device);
    ~WavefrontIntegrator();

    bool RenderSamples(const KernelParams &params);

private:
    static constexpr uint32_t batchSize = 512u * 1024u;

    void AllocateSampleState(WavefrontSampleState &state);
    void FreeSampleState();

    Device *device;
    WavefrontSampleState sampleState{};
    WavefrontSampleState compactionState{};
    std::vector<DeviceBuffer> sampleStateBuffers;
    std::unique_ptr<DeviceKernel> countKernels;
    std::unique_ptr<DeviceKernel> selectSampleIndices;
    std::unique_ptr<DeviceKernel> streamCompaction;
    DeviceBuffer kernelCountsBuffer;
    DeviceBuffer numSampleIndicesBuffer;
    DeviceBuffer numSampleStatesBuffer;
    DeviceBuffer sampleIndicesBuffer;
    DeviceBuffer kernelParamsBuffer;
};

} // namespace ybi
