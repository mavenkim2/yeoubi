#pragma once

#include "device/device_memory_view.h"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>

namespace ybi
{

struct Scene;
enum class RenderKernelId : uint8_t;
struct DispatchParams;
using DevicePtr = uint64_t;
using BVHHandle = uint64_t;

enum class DeviceKind : uint32_t
{
    CPU = 0,
    GPU = 1,
};

enum class DeviceTextureWrapMode : uint8_t
{
    Unknown = 0,
    Repeat = 1,
    Clamp = 2,
    Mirror = 3,
    Black = 4,
    UseMetadata = 5,
};

enum class DeviceTextureFilterMode : uint8_t
{
    Nearest = 0,
    Linear = 1,
};

enum class DeviceTextureFormat : uint8_t
{
    RGBA8_UNORM = 0,
    RGBA16_FLOAT = 1,
};

inline size_t DeviceTextureFormatPixelBytes(DeviceTextureFormat format)
{
    switch (format)
    {
        case DeviceTextureFormat::RGBA8_UNORM:
            return 4u;
        case DeviceTextureFormat::RGBA16_FLOAT:
            return 8u;
        default:
            return 0u;
    }
}

struct DeviceTextureCreateInfo
{
    const void *pixels = nullptr;
    size_t pixelBytes = 0;
    uint32_t width = 0;
    uint32_t height = 0;
    DeviceTextureWrapMode wrapS = DeviceTextureWrapMode::Repeat;
    DeviceTextureWrapMode wrapT = DeviceTextureWrapMode::Repeat;
    DeviceTextureFilterMode filter = DeviceTextureFilterMode::Nearest;
    DeviceTextureFormat format = DeviceTextureFormat::RGBA8_UNORM;
};

struct DeviceTexture
{
    uint64_t handle = 0;
    uint64_t allocation = 0;
    uint32_t width = 0;
    uint32_t height = 0;
    DeviceTextureWrapMode wrapS = DeviceTextureWrapMode::Repeat;
    DeviceTextureWrapMode wrapT = DeviceTextureWrapMode::Repeat;
    DeviceTextureFilterMode filter = DeviceTextureFilterMode::Nearest;
    DeviceTextureFormat format = DeviceTextureFormat::RGBA8_UNORM;
    bool valid = false;
};

struct Device
{
    virtual ~Device() = default;

    static Device *CreateDevice(DeviceKind kind, std::unique_ptr<Device> &storage);

    virtual DeviceKind GetKind() const = 0;
    virtual void BuildBVH(Scene *scene) = 0;
    virtual bool InitializeKernels(const std::string &kernelBlob, std::string *outError) = 0;
    virtual void DestroyKernels() = 0;
    virtual void ClearTransientMemory() = 0;
    virtual bool DispatchKernel(RenderKernelId kernel, const DispatchParams &params) = 0;

    virtual DeviceMemoryView<uint8_t> AllocBytes(size_t numBytes) = 0;
    virtual void FreeBytes(DeviceMemoryView<uint8_t> &view) = 0;
    virtual void CopyBytesToDevice(DeviceMemoryView<uint8_t> dst,
                                   const void *src,
                                   size_t numBytes) = 0;
    virtual void CopyBytesToHost(void *dst,
                                 DeviceMemoryView<const uint8_t> src,
                                 size_t numBytes) = 0;
    virtual bool CreateTexture(const DeviceTextureCreateInfo &info,
                               DeviceTexture *outTexture,
                               std::string *outError) = 0;
    virtual bool UpdateTextureRegion(const DeviceTexture &texture,
                                     uint32_t x,
                                     uint32_t y,
                                     uint32_t width,
                                     uint32_t height,
                                     const void *pixels,
                                     size_t pixelBytes,
                                     std::string *outError) = 0;
    virtual void DestroyTexture(DeviceTexture &texture) = 0;
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

} // namespace ybi
