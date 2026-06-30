#pragma once

#include "device/device.h"
#include "scene/scene.h"
#include <cstdint>
#include <cuda.h>
#include <optix_stubs.h>
#include <string>
#include <vector>

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
    struct ImportedExternalBuffer
    {
        CUdeviceptr ptr = 0;
        size_t size = 0;
        CUexternalMemory externalMemory = nullptr;
    };

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
    std::vector<ImportedExternalBuffer> importedExternalBuffers;

    CUDADevice();
    ~CUDADevice() override;
    DeviceBuffer Alloc(size_t bytes) override;
    void Free(DeviceBuffer &buffer) override;
    DeviceExternalBuffer ImportExternalBuffer(const DisplayBuffer &buffer) override;
    void FreeExternalBuffer(DeviceExternalBuffer &buffer) override;
    DeviceExternalSemaphore ImportExternalSemaphore(const DisplaySemaphore &semaphore) override;
    void FreeExternalSemaphore(DeviceExternalSemaphore &semaphore) override;
    void SignalExternalSemaphore(DeviceExternalSemaphore &semaphore) override;
    BVHHandle BuildTriangleBVH(DeviceArena *deviceArena, Mesh &mesh) override;
    BVHHandle BuildInstancesBVH(DeviceArena *deviceArena,
                                const Array<Instance> &instances,
                                DataContainer *container) override;
    std::unique_ptr<DeviceArena> CreateArena() override;
    void CopyToHost(void *dst, void *ptr, size_t size) override;
    void CopyToHostAsync(void *dst, void *ptr, size_t size) override;
    void CopyToDevice(void *dst, void *ptr, size_t size) override;
    bool DispatchKernel(int width, int height, DeviceBuffer kernelParams) override;

    DeviceTexture CreateTexture(const DeviceTextureInfo &info) override;

    void InitializeOptix();
    bool CreateOptixPipeline(const std::string &ptx);
    void DestroyOptixPipeline();
};

} // namespace ybi
