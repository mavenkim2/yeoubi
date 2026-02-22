#pragma once

#include "device/device_memory_view.h"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>

YBI_NAMESPACE_BEGIN

struct Scene;

enum class DeviceKind : uint32_t
{
    CPU = 0,
    GPU = 1,
};

struct Device
{
    virtual ~Device() = default;

    static Device *CreateDevice(DeviceKind kind, std::unique_ptr<Device> &storage);

    virtual DeviceKind GetKind() const = 0;
    virtual void BuildBVH(Scene *scene) = 0;

    virtual DeviceMemoryView<uint8_t> AllocBytes(size_t numBytes) = 0;
    virtual void FreeBytes(DeviceMemoryView<uint8_t> &view) = 0;
    virtual void CopyBytesToDevice(DeviceMemoryView<uint8_t> dst,
                                   const void *src,
                                   size_t numBytes) = 0;
    virtual void CopyBytesToHost(void *dst,
                                 DeviceMemoryView<const uint8_t> src,
                                 size_t numBytes) = 0;
    virtual size_t GetBVHAllocatedBytes() const = 0;

    template <typename T>
    DeviceMemoryView<T> Alloc(size_t count)
    {
        if (count == 0)
        {
            return {};
        }
        DeviceMemoryView<uint8_t> bytes = AllocBytes(sizeof(T) * count);
        return {(T *)bytes.data(), count};
    }

    template <typename T>
    void Free(DeviceMemoryView<T> &view)
    {
        if (view.data() == nullptr || view.size() == 0)
        {
            return;
        }
        DeviceMemoryView<uint8_t> bytes = ByteView(view);
        FreeBytes(bytes);
        view = {};
    }

    template <typename T>
    static DeviceMemoryView<uint8_t> ByteView(DeviceMemoryView<T> view)
    {
        return {(uint8_t *)view.data(), view.numBytes()};
    }

    template <typename T>
    static DeviceMemoryView<const uint8_t> ByteView(DeviceMemoryView<const T> view)
    {
        return {(const uint8_t *)view.data(), view.numBytes()};
    }
};

YBI_NAMESPACE_END
