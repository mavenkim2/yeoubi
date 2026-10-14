#pragma once

#include "device/device.h"
#include "kernel/guiding/types.h"

namespace ybi
{

class PathGuidingSystem
{
    Device *device;
    DeviceBuffer kernelParams;
    DeviceBuffer activeNodeItems;
    uint32_t generation = 0;

    std::unique_ptr<DeviceKernel> generateWorkItems;
    std::unique_ptr<DeviceKernel> initializeNodePositionStatistics;
    std::unique_ptr<DeviceKernel> calculateNodeBounds;
    std::unique_ptr<DeviceKernel> calculateNodePositionStatistics;
    std::unique_ptr<DeviceKernel> selectNodeSplits;
    std::unique_ptr<DeviceKernel> partitionCount;
    std::unique_ptr<DeviceKernel> partitionExclusiveScan;
    std::unique_ptr<DeviceKernel> partitionSamples;

public:
    PathGuidingSystem(Device *device, DeviceKind deviceKind);
    ~PathGuidingSystem();
    PathGuidingSystem(const PathGuidingSystem &) = delete;
    PathGuidingSystem &operator=(const PathGuidingSystem &) = delete;

    void UpdateKDTree(PathGuidingKdTreeParams &params);
};

} // namespace ybi
