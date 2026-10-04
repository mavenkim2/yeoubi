#pragma once

#include "device/device.h"

namespace ybi 
{

struct CPUDevice : Device
{
    size_t totalAllocated;
    size_t bvhTotalAllocated;
    bool supportsClusterBVHs;

    explicit CPUDevice(bool createRendererPipeline = true);
    ~CPUDevice() override;
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
    bool DispatchKernel(KernelType kernelType,
                        uint32_t numWorkItems,
                        uint32_t blockSize,
                        DeviceBuffer kernelParams,
                        float *elapsedMilliseconds = nullptr) override;

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
    void ExecuteKernel(DeviceKernel &kernel,
                       const DeviceKernelLaunchDesc &desc,
                       float *elapsedMilliseconds = nullptr) override;

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
};

}
