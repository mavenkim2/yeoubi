#pragma once

#include "device/device.h"
#include "kernel/guiding/types.h"

namespace ybi
{

class PathGuidingSystem
{
    Device *device;
    DeviceBuffer kernelParams;

    std::unique_ptr<DeviceKernel> generateWorkItems;
    std::unique_ptr<DeviceKernel> initializeNodePositionStatistics;
    std::unique_ptr<DeviceKernel> calculateNodeBounds;
    std::unique_ptr<DeviceKernel> calculateNodePositionStatistics;

public:
    PathGuidingSystem(Device *device, DeviceKind deviceKind);
    ~PathGuidingSystem();
    PathGuidingSystem(const PathGuidingSystem &) = delete;
    PathGuidingSystem &operator=(const PathGuidingSystem &) = delete;

    void UpdateKDTree(PathGuidingKdTreeParams &params);
};

} // namespace ybi
