#pragma once

#include "core/vec2u.h"
#include "device/device_memory.h"
#include "kernel/state/types.h"
#include <memory>
#include <string>
#include <vector>

namespace ybi
{

struct Device;
struct DeviceKernel;
struct KernelParams;
struct WavefrontParams;
#if defined(YBI_ENABLE_PATH_GUIDING)
class PathGuidingSystem;
#endif

class WavefrontIntegrator
{
public:
    WavefrontIntegrator(Device *device);
    ~WavefrontIntegrator();

    bool RenderSamples(const KernelParams &params);
    bool RenderSamples(const KernelParams &params, Vec2u outputSize);
    bool ExportGuidingCells(const std::string &path, std::string *error = nullptr) const;
    void FinishGuidingTraining();

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
#if defined(YBI_ENABLE_PATH_GUIDING)
    std::unique_ptr<PathGuidingSystem> pathGuidingSystem;
#endif
};

} // namespace ybi
