#pragma once

#include "device/device.h"
#include "scene/scene.h"
#include <cstdint>
#include <optix_stubs.h>
#include <string>

namespace ybi
{

struct DeviceBVH
{
    BVHHandle handle;
};

struct CUDAOptixPipeline
{
    OptixModule module = nullptr;
    OptixProgramGroup raygenProgramGroup = nullptr;
    OptixProgramGroup missProgramGroup = nullptr;
    OptixProgramGroup hitgroupProgramGroup = nullptr;
    OptixPipeline pipeline = nullptr;
    CUdeviceptr raygenRecord = 0;
    CUdeviceptr missRecord = 0;
    CUdeviceptr hitgroupRecord = 0;
    OptixShaderBindingTable shaderBindingTable = {};
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
    CUDAOptixPipeline optixPipeline;

    CUDADevice();
    ~CUDADevice() override;
    DeviceBuffer Alloc(size_t bytes) override;
    void Free(DeviceBuffer &buffer) override;
    BVHHandle BuildTriangleBVH(DeviceArena *deviceArena, Mesh &mesh) override;
    BVHHandle BuildInstancesBVH(DeviceArena *deviceArena,
                                const Array<Instance> &instances,
                                DataContainer *container) override;
    std::unique_ptr<DeviceArena> CreateArena() override;
    void CopyToDevice(void *dst, void *ptr, size_t size) override;
    bool DispatchKernel(int width, int height, DeviceBuffer kernelParams) override;

    void InitializeOptix();
    bool CreateOptixPipeline(const std::string &ptx);
    void DestroyOptixPipeline();
};

} // namespace ybi
