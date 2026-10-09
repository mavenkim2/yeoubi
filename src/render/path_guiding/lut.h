#pragma once

#include "device/device_texture.h"
#include "kernel/guiding/types.h"

namespace ybi
{

struct Device;

constexpr const char *VMFSplitTableFilename = "vmf_split_table.bin";

// File: float maxInvSqrtK, then 256 float cos(beta) values.
struct VMFSplitTable
{
    DeviceTexture texture{};
    float maxInvSqrtK = 0.f;

    void Bind(GuidingParams &params) const
    {
        params.vmfSplitTexture = texture.textureObject;
        params.vmfSplitMaxInvSqrtK = maxInvSqrtK;
    }
};

VMFSplitTable LoadVMFSplitTable(Device &device);

} // namespace ybi
