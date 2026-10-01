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

struct CUDADeviceKernel : DeviceKernel
{
    CUcontext context;
    CUmodule module;
    CUfunction function;

    CUDADeviceKernel(CUcontext context, CUmodule module, CUfunction function);
    ~CUDADeviceKernel() override;
};

struct CUDADeviceStream : DeviceStream
{
    CUcontext context = nullptr;
    CUstream stream = nullptr;

    explicit CUDADeviceStream(CUcontext context);
    ~CUDADeviceStream() override;
};

struct CUDADeviceEvent : DeviceEvent
{
    CUcontext context = nullptr;
    CUevent event = nullptr;

    explicit CUDADeviceEvent(CUcontext context);
    ~CUDADeviceEvent() override;
};

struct CUDAOptixPipeline
{
    OptixModule module = nullptr;
    OptixProgramGroup raygenProgramGroup = nullptr;
    OptixProgramGroup missProgramGroup = nullptr;
    OptixProgramGroup hitgroupProgramGroup = nullptr;
    std::vector<OptixModule> callableModules;
    std::vector<OptixProgramGroup> callableProgramGroups;
    OptixPipeline pipeline = nullptr;
    CUdeviceptr raygenRecord = 0;
    CUdeviceptr missRecord = 0;
    CUdeviceptr hitgroupRecord = 0;
    CUdeviceptr callableRecords = 0;
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
    bool supportsClusterBVHs;

    explicit CUDADevice(bool createRendererPipeline = true);
    ~CUDADevice() override;
    void WaitIdle() override;
    DeviceBuffer Alloc(size_t bytes) override;
    void Free(DeviceBuffer &buffer) override;
    DeviceBuffer AllocHost(size_t bytes) override;
    void FreeHost(DeviceBuffer &buffer) override;
    DeviceExternalBuffer ImportExternalBuffer(const DisplayBuffer &buffer) override;
    void FreeExternalBuffer(DeviceExternalBuffer &buffer) override;
    DeviceExternalSemaphore ImportExternalSemaphore(const DisplaySemaphore &semaphore) override;
    void FreeExternalSemaphore(DeviceExternalSemaphore &semaphore) override;
    void SignalExternalSemaphore(DeviceExternalSemaphore &semaphore) override;
    BVHHandle
    BuildTriangleBVH(DeviceArena *deviceArena, DataContainer *container, Mesh &mesh) override;
    BVHHandle BuildInstancesBVH(DeviceArena *deviceArena,
                                const Array<Instance> &instances,
                                DataContainer *container,
                                const Array<int> &objectToBVHIndex) override;
    std::unique_ptr<DeviceArena> CreateArena() override;
    void CopyToHost(void *dst, void *ptr, size_t size) override;
    void CopyToHostAsync(void *dst, void *ptr, size_t size) override;
    void CopyToDevice(void *dst, void *ptr, size_t size) override;
    void CopyToDeviceAsync(void *dst, void *ptr, size_t size) override;
    void Memzero(void *dst, size_t size) override;
    bool DispatchKernel(int width, int height, DeviceBuffer kernelParams) override;

    DeviceTexture CreateTexture(const DeviceTextureInfo &info) override;
    void DestroyTexture(DeviceTexture &texture) override;
    std::unique_ptr<DeviceStream> CreateDeviceStream() override;
    std::shared_ptr<DeviceEvent> CreateDeviceEvent() override;
    void SyncDeviceStream(DeviceStream &stream) override;
    void RecordDeviceEvent(DeviceStream *stream, DeviceEvent &event) override;
    void WaitDeviceEvent(DeviceStream &stream, DeviceEvent &event) override;
    void CopyToDevice3DAsync(DeviceStream &stream,
                             DeviceTexture &texture,
                             const DeviceTextureCopyDesc &desc) override;
    void CopyToDevice3DBatchAsync(DeviceStream &stream,
                                  const DeviceTextureBatchCopyDesc *descs,
                                  size_t numDescs) override;
    std::unique_ptr<DeviceKernel> CreateDeviceKernel(const DeviceKernelDesc &desc) override;
    bool CompileKernel(const std::string &source, DeviceModule &deviceModule) override;
    bool RegisterDirectCallables(DeviceModule deviceModule,
                                 const std::vector<std::string> &entryNames,
                                 std::vector<uint32_t> &callableIndices) override;
    std::string GetCallableFunctionSignature() const override;
    std::string GetDirectCallablePrefix() const override;
    void ExecuteKernel(DeviceKernel &kernel, const DeviceKernelLaunchDesc &desc) override;

    void InitializeOptix(bool createRendererPipeline);
    bool CreateOptixPipeline(const std::string &ptx);
    void DestroyOptixPipeline();

    // CLAS BVHs
    bool SupportsClusterBVHs() const override;
    BVHHandle BuildClusterBVH(DeviceArena *deviceArena,
                              DataContainer *container,
                              const Mesh &mesh,
                              const MeshletClustersResult &clusterResult) override;
    size_t GetBVHAllocationSize() const override;
    void StartProfiler() const override;
    void StopProfiler() const override;
    void CreatePipeline() override;
    void DestroyPipeline() override;
    bool CreatePipeline(MaterialBackend backend) override;

private:
    void GetOptixCompileOptions(OptixModuleCompileOptions &moduleOptions,
                                OptixPipelineCompileOptions &pipelineOptions);
    bool LinkOptixPipeline(const OptixPipelineCompileOptions &pipelineCompileOptions);
    void DestroyLinkedOptixPipeline();

};

} // namespace ybi
