#pragma once

#include "device/device_memory.h"
#include <cstddef>
#include <cstdint>
#include <memory>

struct GLFWwindow;

namespace ybi
{

enum class ExternalMemoryHandleType
{
    None,
    OpaqueWin32,
    OpaqueFd,
    D3D12Resource,
    D3D12Heap,
};

struct ExternalMemoryHandle
{
    ExternalMemoryHandleType type = ExternalMemoryHandleType::None;
    void *handle = nullptr;
    int fd = -1;
    size_t size = 0;
};

struct DisplayBuffer
{
    virtual ~DisplayBuffer() = default;

    size_t size = 0;
    ExternalMemoryHandle memory;
};

struct DeviceExternalBuffer
{
    void *handle = nullptr;
    DeviceBuffer buffer;
};

struct DisplaySemaphore
{
    virtual ~DisplaySemaphore() = default;

    ExternalMemoryHandle handle;
};

struct DeviceExternalSemaphore
{
    void *handle = nullptr;
};

struct DisplayDevice
{
    uint32_t numFramesInFlight = 2;
    uint32_t currentFrame = 0;
    uint64_t frameCount = 0;

    virtual ~DisplayDevice() = default;

    static std::unique_ptr<DisplayDevice> CreateDevice(GLFWwindow *window);
    virtual std::unique_ptr<DisplayBuffer> AllocateExternalBuffer(size_t size) = 0;
    virtual std::unique_ptr<DisplaySemaphore> CreateDisplaySemaphore() = 0;
    // Complete submitted display work before releasing resources it references.
    virtual void WaitIdle() = 0;
    virtual bool BeginFrame() = 0;
    virtual bool
    EndFrame(DisplayBuffer &framebuffer, DisplaySemaphore &cudaFinishedSemaphore, int width, int height) = 0;

    uint32_t GetNumFramesInFlight() const
    {
        return numFramesInFlight;
    }

    uint32_t GetCurrentFrameIndex() const
    {
        return currentFrame;
    }
};

} // namespace ybi
