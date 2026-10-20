#pragma once

#include "core/array.h"
#include "core/assert.h"
#include "device/device_memory.h"
#include "device/device_texture.h"
#include "display_device/display_device.h"
#include "kernel/state/types.h"
#include "scene/bvh.h"
#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <type_traits>
#include <vector>

namespace ybi
{

class Arena;
struct DataContainer;
struct DisplayBuffer;
struct Instance;
class Mesh;
struct MeshletClustersResult;

enum class DeviceKind
{
    CPU,
    GPU,
};

enum class MaterialBackend
{
    OpenPbr,
    MaterialX,
};

class DeviceArena
{
protected:
    struct Mark
    {
        DeviceArena *arena;
        DeviceArena *current;
        size_t used;
    };

    virtual void RestoreMark(Mark mark) = 0;
    virtual DeviceArena *GetCurrent() const = 0;
    virtual size_t GetUsed() const = 0;

public:
    virtual ~DeviceArena() = default;

    struct Scope
    {
        Mark mark;

        explicit Scope(DeviceArena *arena)
        {
            YBI_ASSERT(arena);
            mark.arena = arena;
            mark.current = arena->GetCurrent();
            mark.used = mark.current->GetUsed();
        }
        ~Scope()
        {
            if (mark.arena)
            {
                mark.arena->RestoreMark(mark);
            }
        }
        Scope(const Scope &) = delete;
        Scope &operator=(const Scope &) = delete;
    };
};

typedef uint64_t DeviceModule;
struct DeviceKernel
{
    virtual ~DeviceKernel() = default;
};

struct DeviceStream
{
    virtual ~DeviceStream() = default;
};

struct DeviceEvent
{
    virtual ~DeviceEvent() = default;
};

struct DeviceKernelDesc
{
    std::string filename;
    std::string entryName;
};

struct DeviceKernelLaunchDesc
{
    uint32_t gridX = 1;
    uint32_t gridY = 1;
    uint32_t gridZ = 1;
    uint32_t blockX = 1;
    uint32_t blockY = 1;
    uint32_t blockZ = 1;
    uint32_t sharedMemoryBytes = 0;
    void **arguments = nullptr;
};

struct DeviceTextureCopyDesc
{
    const void *src = nullptr;
    size_t srcRowPitch = 0;
    size_t srcSlicePitch = 0;
    size_t widthInBytes = 0;
    size_t height = 0;
    size_t depth = 1;
    // Destination coordinates are in texels. Backends convert to API-specific units.
    uint32_t dstX = 0;
    uint32_t dstY = 0;
    uint32_t dstZ = 0;
};

struct DeviceTextureBatchCopyDesc
{
    DeviceTexture *texture = nullptr;
    DeviceTextureCopyDesc copy;
};

struct Device
{
    virtual ~Device() = default;

    virtual void WaitIdle() = 0;

    virtual DeviceBuffer<std::byte> AllocBytes(size_t bytes) = 0;
    virtual void FreeBytes(DeviceBuffer<std::byte> &buffer) = 0;
    virtual DeviceBuffer<std::byte> AllocHostBytes(size_t bytes) = 0;
    virtual void FreeHostBytes(DeviceBuffer<std::byte> &buffer) = 0;

    template <typename T = std::byte>
    DeviceBuffer<T> Alloc(size_t count)
    {
        return AllocBytes(count * sizeof(T)).template As<T>();
    }

    template <typename T>
    void Free(DeviceBuffer<T> &buffer)
    {
        auto bytes = buffer.template As<std::byte>();
        FreeBytes(bytes);
        buffer = {};
    }

    template <typename T = std::byte>
    DeviceBuffer<T> AllocHost(size_t count)
    {
        return AllocHostBytes(count * sizeof(T)).template As<T>();
    }

    template <typename T>
    void FreeHost(DeviceBuffer<T> &buffer)
    {
        auto bytes = buffer.template As<std::byte>();
        FreeHostBytes(bytes);
        buffer = {};
    }

    virtual DeviceExternalBuffer ImportExternalBuffer(const DisplayBuffer &buffer) = 0;
    virtual void FreeExternalBuffer(DeviceExternalBuffer &buffer) = 0;
    virtual DeviceExternalSemaphore ImportExternalSemaphore(const DisplaySemaphore &semaphore) = 0;
    virtual void FreeExternalSemaphore(DeviceExternalSemaphore &semaphore) = 0;
    virtual void SignalExternalSemaphore(DeviceExternalSemaphore &semaphore) = 0;
    virtual BVHHandle
    BuildTriangleBVH(DeviceArena *deviceArena, DataContainer *container, Mesh &mesh) = 0;
    virtual BVHHandle BuildInstancesBVH(DeviceArena *deviceArena,
                                        const Array<Instance> &instances,
                                        DataContainer *container,
                                        const Array<int> &objectToBVHIndex) = 0;
    virtual std::unique_ptr<DeviceArena> CreateArena() = 0;
    static std::unique_ptr<Device> CreateDevice(DeviceKind kind,
                                                bool createRendererPipeline = true);
    virtual void CopyToDevice(void *dst, void *ptr, size_t size) = 0;
    virtual void CopyToDeviceAsync(void *dst, void *ptr, size_t size) = 0;
    virtual void CopyToHost(void *dst, void *ptr, size_t size) = 0;
    virtual void CopyToHostAsync(void *dst, void *ptr, size_t size) = 0;
    virtual void Memzero(void *dst, size_t size) = 0;
    virtual void MemzeroAsync(void *dst, size_t size) = 0;
    // Enqueues work on the default stream; synchronize when results are needed.
    virtual bool DispatchKernel(KernelType kernelType,
                                uint32_t numWorkItems,
                                uint32_t blockSize,
                                DeviceBuffer<std::byte> kernelParams) = 0;

    template <typename T>
    bool DispatchKernel(KernelType kernelType,
                        uint32_t numWorkItems,
                        uint32_t blockSize,
                        DeviceBuffer<T> kernelParams)
    {
        return DispatchKernel(
            kernelType, numWorkItems, blockSize, kernelParams.template As<std::byte>());
    }

    virtual DeviceTexture CreateTexture(const DeviceTextureInfo &info) = 0;
    virtual void DestroyTexture(DeviceTexture &texture) = 0;
    virtual std::unique_ptr<DeviceStream> CreateDeviceStream() = 0;
    virtual std::shared_ptr<DeviceEvent> CreateDeviceEvent() = 0;
    virtual void SyncDeviceStream(DeviceStream *stream = nullptr) = 0;
    void SyncDeviceStream(DeviceStream &stream) { SyncDeviceStream(&stream); }
    virtual void RecordDeviceEvent(DeviceStream *stream, DeviceEvent &event) = 0;
    virtual void WaitDeviceEvent(DeviceStream &stream, DeviceEvent &event) = 0;
    virtual void CopyToDevice3D(DeviceTexture &texture, const DeviceTextureCopyDesc &desc) = 0;
    virtual void CopyToDevice3DAsync(DeviceStream &stream,
                                     DeviceTexture &texture,
                                     const DeviceTextureCopyDesc &desc) = 0;
    virtual void CopyToDevice3DBatchAsync(DeviceStream &stream,
                                          const DeviceTextureBatchCopyDesc *descs,
                                          size_t numDescs) = 0;
    virtual std::unique_ptr<DeviceKernel> CreateDeviceKernel(const DeviceKernelDesc &desc) = 0;
    virtual bool CompileKernel(const std::string &source, DeviceModule &deviceModule) = 0;
    virtual bool RegisterDirectCallables(DeviceModule deviceModule,
                                         const std::vector<std::string> &entryNames,
                                         std::vector<uint32_t> &callableIndices) = 0;
    virtual std::string GetCallableFunctionSignature() const = 0;
    virtual std::string GetDirectCallablePrefix() const = 0;
    virtual void ExecuteKernel(DeviceKernel &kernel,
                               const DeviceKernelLaunchDesc &desc) = 0;

    virtual bool SupportsClusterBVHs() const = 0;
    virtual BVHHandle BuildClusterBVH(DeviceArena *deviceArena,
                                      DataContainer *container,
                                      const Mesh &mesh,
                                      const MeshletClustersResult &clusterResult) = 0;

    virtual size_t GetBVHAllocationSize() const = 0;
    virtual void StartProfiler() const = 0;
    virtual void StopProfiler() const = 0;

    virtual void CreatePipeline() = 0;
    virtual void DestroyPipeline() {}
    virtual bool CreatePipeline(MaterialBackend backend)
    {
        if (backend != MaterialBackend::OpenPbr)
            return false;
        CreatePipeline();
        return true;
    }
    virtual void CopyDeviceToDeviceAsync(void *dst, void *ptr, size_t size) = 0;
};

} // namespace ybi
