#include "display_device/display_device.h"
#include "device/device.h"
#include <utility>

namespace ybi
{

DeviceExternalBuffer::DeviceExternalBuffer(DeviceExternalBuffer &&other) noexcept
    : handle(std::exchange(other.handle, nullptr)), device(std::exchange(other.device, nullptr)),
      buffer(std::exchange(other.buffer, DeviceBuffer{}))
{
}

DeviceExternalBuffer &DeviceExternalBuffer::operator=(DeviceExternalBuffer &&other)
{
    if (this != &other)
    {
        if (device)
            device->FreeExternalBuffer(*this);
        handle = std::exchange(other.handle, nullptr);
        device = std::exchange(other.device, nullptr);
        buffer = std::exchange(other.buffer, DeviceBuffer{});
    }
    return *this;
}

DeviceExternalBuffer::~DeviceExternalBuffer()
{
    if (device)
        device->FreeExternalBuffer(*this);
}

DeviceExternalSemaphore::DeviceExternalSemaphore(DeviceExternalSemaphore &&other) noexcept
    : handle(std::exchange(other.handle, nullptr)), device(std::exchange(other.device, nullptr))
{
}

DeviceExternalSemaphore &DeviceExternalSemaphore::operator=(DeviceExternalSemaphore &&other)
{
    if (this != &other)
    {
        if (device)
            device->FreeExternalSemaphore(*this);
        handle = std::exchange(other.handle, nullptr);
        device = std::exchange(other.device, nullptr);
    }
    return *this;
}

DeviceExternalSemaphore::~DeviceExternalSemaphore()
{
    if (device)
        device->FreeExternalSemaphore(*this);
}

} // namespace ybi
