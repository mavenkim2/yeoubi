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
    DeviceBuffer nodes;
    DeviceBuffer nodeItems;
    DeviceBuffer nextNodeItems;
    DeviceBuffer workItemLeftOffsets;
    uint32_t generation = 0;

    std::unique_ptr<DeviceKernel> generateNodeItems;
    std::unique_ptr<DeviceKernel> generateWorkItems;
    std::unique_ptr<DeviceKernel> initializeNodePositionStatistics;
    std::unique_ptr<DeviceKernel> calculateNodeBounds;
    std::unique_ptr<DeviceKernel> calculateNodePositionStatistics;
    std::unique_ptr<DeviceKernel> selectNodeSplits;
    std::unique_ptr<DeviceKernel> partitionCount;
    std::unique_ptr<DeviceKernel> partitionExclusiveScan;
    std::unique_ptr<DeviceKernel> partitionSamples;

    void ReserveBuffer(DeviceBuffer &buffer, size_t size, size_t preserveSize = 0);
    void CopyToOwnedBuffer(DeviceBuffer &buffer, void *source, size_t size);
    bool PartitionSamples(PathGuidingKdTreeParams &params);

public:
    PathGuidingSystem(Device *device, DeviceKind deviceKind);
    ~PathGuidingSystem();
    PathGuidingSystem(const PathGuidingSystem &) = delete;
    PathGuidingSystem &operator=(const PathGuidingSystem &) = delete;

    void UpdateKDTree(PathGuidingKdTreeParams &params);
};

} // namespace ybi
