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
struct WavefrontParams;
#if defined(YBI_USE_PATH_GUIDING)
struct PathGuidingTrainingSample;
#endif

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
    std::vector<DeviceBuffer<std::byte>> sampleStateBuffers;
    std::unique_ptr<DeviceKernel> countKernels;
    std::unique_ptr<DeviceKernel> selectSampleIndices;
    std::unique_ptr<DeviceKernel> streamCompaction;
    DeviceBuffer<KernelCounts> kernelCountsBuffer;
    DeviceBuffer<KernelCounts> hostKernelCountsBuffer;
    DeviceBuffer<std::byte> hostUploadBuffer;
    DeviceBuffer<uint32_t> numSampleIndicesBuffer;
    DeviceBuffer<uint32_t> numSampleStatesBuffer;
    DeviceBuffer<uint32_t> sampleIndicesBuffer;
    DeviceBuffer<WavefrontParams> kernelParamsBuffer;
#if defined(YBI_USE_PATH_GUIDING)
    DeviceBuffer<uint32_t> radianceRecordAllocationCursorBuffer;
    DeviceBuffer<uint32_t> radianceRecordRetiredCountBuffer;
    DeviceBuffer<PathGuidingTrainingSample> pathGuidingTrainingSamplesBuffer;
    DeviceBuffer<uint32_t> pathGuidingTrainingSampleCountBuffer;
#endif
};

} // namespace ybi
