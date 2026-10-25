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
    const char *name = nullptr;

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

    CUDADeviceEvent(CUcontext context, bool enableTiming = false);
    ~CUDADeviceEvent() override;
};

enum OptixProgramGroups
{
    PROGRAM_GROUP_RAYGEN_INTERSECT,
    PROGRAM_GROUP_RAYGEN_SHADE_HIT,
    PROGRAM_GROUP_RAYGEN_DIRECT_LIGHTING,
    PROGRAM_GROUP_RAYGEN_INDIRECT_LIGHTING,
    PROGRAM_GROUP_RAYGEN_END,

    PROGRAM_GROUP_MISS = PROGRAM_GROUP_RAYGEN_END,
    PROGRAM_GROUP_CLOSEST_HIT,
    PROGRAM_GROUP_COUNT,
};

struct alignas(OPTIX_SBT_RECORD_ALIGNMENT) PipelineSbtRecord
{
    char header[OPTIX_SBT_RECORD_HEADER_SIZE] = {};
};

struct CUDAOptixPipeline
{
    OptixModule module = nullptr;
    OptixProgramGroup programGroups[PROGRAM_GROUP_COUNT] = {};
    PipelineSbtRecord sbtRecords[PROGRAM_GROUP_COUNT] = {};
    CUdeviceptr sbtRecordsBuffer = 0;
    std::vector<OptixModule> callableModules;
    std::vector<OptixProgramGroup> callableProgramGroups;
    OptixPipeline pipeline = nullptr;
    CUdeviceptr callableRecords = 0;
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
    DeviceBuffer<std::byte> AllocBytes(size_t bytes) override;
    void FreeBytes(DeviceBuffer<std::byte> &buffer) override;
    DeviceBuffer<std::byte> AllocHostBytes(size_t bytes) override;
    void FreeHostBytes(DeviceBuffer<std::byte> &buffer) override;
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
    void CopyDeviceToDeviceAsync(void *dst, void *ptr, size_t size) override;
    void Memzero(void *dst, size_t size) override;
    void MemzeroAsync(void *dst, size_t size) override;
    using Device::DispatchKernel;
    bool DispatchKernel(KernelType kernelType,
                        uint32_t numWorkItems,
                        uint32_t blockSize,
                        DeviceBuffer<std::byte> kernelParams) override;

    DeviceTexture CreateTexture(const DeviceTextureInfo &info) override;
    void DestroyTexture(DeviceTexture &texture) override;
    std::unique_ptr<DeviceStream> CreateDeviceStream() override;
    std::shared_ptr<DeviceEvent> CreateDeviceEvent() override;
    std::unique_ptr<DeviceEventTimer> StartEventTimer(DeviceStream *stream = nullptr) override;
    void StopEventTimer(DeviceEventTimer &timer, DeviceStream *stream = nullptr) override;
    using Device::SyncDeviceStream;
    void SyncDeviceStream(DeviceStream *stream = nullptr) override;
    void RecordDeviceEvent(DeviceStream *stream, DeviceEvent &event) override;
    void WaitDeviceEvent(DeviceStream &stream, DeviceEvent &event) override;
    void CopyToDevice3D(DeviceTexture &texture, const DeviceTextureCopyDesc &desc) override;
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
    void ExecuteKernel(DeviceKernel &kernel,
                       const DeviceKernelLaunchDesc &desc) override;

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
    std::unique_ptr<DeviceKernel> generateCameraRays;
    void GetOptixCompileOptions(OptixModuleCompileOptions &moduleOptions,
                                OptixPipelineCompileOptions &pipelineOptions);
    bool LinkOptixPipeline(const OptixPipelineCompileOptions &pipelineCompileOptions);
    void DestroyLinkedOptixPipeline();
};

} // namespace ybi
