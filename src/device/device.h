#pragma once

#include "core/array.h"
#include "core/assert.h"
#include "device/device_memory.h"
#include "device/device_texture.h"
#include "display_device/display_device.h"
#include "scene/bvh.h"
#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <type_traits>

namespace ybi
{

class Arena;
struct DataContainer;
struct DisplayBuffer;
struct Instance;
class Mesh;

enum class DeviceKind
{
    CPU,
    GPU,
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
    size_t widthInBytes = 0;
    size_t height = 0;
    // Destination coordinates are in texels. Backends convert to API-specific units.
    uint32_t dstX = 0;
    uint32_t dstY = 0;
    uint32_t dstLayer = 0;
};

struct DeviceTextureBatchCopyDesc
{
    DeviceTexture *texture = nullptr;
    DeviceTextureCopyDesc copy;
};

struct Device
{
    virtual ~Device() = default;

    virtual DeviceBuffer Alloc(size_t bytes) = 0;
    virtual void Free(DeviceBuffer &buffer) = 0;
    virtual DeviceBuffer AllocHost(size_t bytes) = 0;
    virtual void FreeHost(DeviceBuffer &buffer) = 0;
    virtual DeviceExternalBuffer ImportExternalBuffer(const DisplayBuffer &buffer) = 0;
    virtual void FreeExternalBuffer(DeviceExternalBuffer &buffer) = 0;
    virtual DeviceExternalSemaphore ImportExternalSemaphore(const DisplaySemaphore &semaphore) = 0;
    virtual void FreeExternalSemaphore(DeviceExternalSemaphore &semaphore) = 0;
    virtual void SignalExternalSemaphore(DeviceExternalSemaphore &semaphore) = 0;
    virtual BVHHandle BuildTriangleBVH(DeviceArena *deviceArena, Mesh &mesh) = 0;
    virtual BVHHandle BuildInstancesBVH(DeviceArena *deviceArena,
                                        const Array<Instance> &instances,
                                        DataContainer *container) = 0;
    virtual std::unique_ptr<DeviceArena> CreateArena() = 0;
    static std::unique_ptr<Device> CreateDevice(DeviceKind kind);
    virtual void CopyToDevice(void *dst, void *ptr, size_t size) = 0;
    virtual void CopyToDeviceAsync(void *dst, void *ptr, size_t size) = 0;
    virtual void CopyToHost(void *dst, void *ptr, size_t size) = 0;
    virtual void CopyToHostAsync(void *dst, void *ptr, size_t size) = 0;
    virtual void Memzero(void *dst, size_t size) = 0;
    virtual bool DispatchKernel(int width, int height, DeviceBuffer kernelParams) = 0;
    virtual DeviceTexture CreateTexture(const DeviceTextureInfo &info) = 0;
    virtual std::unique_ptr<DeviceStream> CreateDeviceStream() = 0;
    virtual std::shared_ptr<DeviceEvent> CreateDeviceEvent() = 0;
    virtual void SyncDeviceStream(DeviceStream &stream) = 0;
    virtual void RecordDeviceEvent(DeviceStream *stream, DeviceEvent &event) = 0;
    virtual void WaitDeviceEvent(DeviceStream &stream, DeviceEvent &event) = 0;
    virtual void CopyToDevice3DAsync(DeviceStream &stream,
                                     DeviceTexture &texture,
                                     const DeviceTextureCopyDesc &desc) = 0;
    virtual void CopyToDevice3DBatchAsync(DeviceStream &stream,
                                          const DeviceTextureBatchCopyDesc *descs,
                                          size_t numDescs) = 0;
    virtual std::unique_ptr<DeviceKernel> CreateDeviceKernel(const DeviceKernelDesc &desc) = 0;
    virtual void ExecuteKernel(DeviceKernel &kernel, const DeviceKernelLaunchDesc &desc) = 0;
};

} // namespace ybi
