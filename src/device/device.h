#pragma once

#include "core/array.h"
#include "core/assert.h"
#include "device/device_memory.h"
#include "scene/bvh.h"
#include <cstddef>
#include <cstdint>
#include <memory>
#include <type_traits>

namespace ybi
{

class Arena;
struct DataContainer;
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

struct Device
{
    virtual ~Device() = default;

    virtual DeviceBuffer Alloc(size_t bytes) = 0;
    virtual void Free(DeviceBuffer &buffer) = 0;
    virtual BVHHandle BuildTriangleBVH(DeviceArena *deviceArena, Mesh &mesh) = 0;
    virtual BVHHandle BuildInstancesBVH(DeviceArena *deviceArena,
                                        const Array<Instance> &instances,
                                        DataContainer *container) = 0;
    virtual std::unique_ptr<DeviceArena> CreateArena() = 0;
    static std::unique_ptr<Device> CreateDevice(DeviceKind kind);
    virtual void CopyToDevice(void *dst, void *ptr, size_t size) = 0;
    virtual void CopyToHost(void *dst, void *ptr, size_t size) = 0;
    virtual bool DispatchKernel(int width, int height, DeviceBuffer kernelParams) = 0;
};

} // namespace ybi
