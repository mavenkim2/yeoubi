#pragma once

#include "device/device.h"
#include "scene/scene.h"
#include <cstdint>
#include <optix_stubs.h>

namespace ybi
{

struct DeviceBVH
{
    BVHHandle handle;
};

struct CUDADevice : Device
{
    struct ContextScope
    {
        ContextScope(CUDADevice *device);
        ~ContextScope();
    };

    CUdevice device;
    CUcontext cudaContext;
    OptixDeviceContext optixDeviceContext;
    size_t totalAllocated;
    size_t bvhTotalAllocated;

    CUDADevice();
    DeviceBuffer Alloc(size_t bytes) override;
    void Free(DeviceBuffer &buffer) override;
    BVHHandle BuildTriangleBVH(DeviceArena *deviceArena, Mesh &mesh) override;
    BVHHandle BuildInstancesBVH(DeviceArena *deviceArena,
                                const Array<Instance> &instances,
                                DataContainer *container) override;
    std::unique_ptr<DeviceArena> CreateArena() override;

    void InitializeOptix();
};

} // namespace ybi
