#pragma once

#include "device/device.h"
#include "kernel/guiding/kd_tree.h"
#include "kernel/guiding/vmf.h"

namespace ybi
{

struct Device;

class PathGuidingSystem
{
    Device *device;
    PathGuidingKdTreeNode *nodes;
    PathGuidingKdTreeNodeItem *nodeItems;

    std::unique_ptr<DeviceKernel> generateWorkItems;
    std::unique_ptr<DeviceKernel> calculateNodeBounds;
    std::unique_ptr<DeviceKernel> calculateNodePositionStatistics;

    PathGuidingSystem(Device *device);
    void Partition();

public:
    void UpdateKDTree();
};

} // namespace ybi
