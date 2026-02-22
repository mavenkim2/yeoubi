#include "device/cpu_device.h"

#if defined(WITH_EMBREE)

#include <cstdio>
#include <cstdlib>
#include <cstring>

YBI_NAMESPACE_BEGIN

CPUDevice::CPUDevice()
{
    embreeDevice = rtcNewDevice(nullptr);
    if (!embreeDevice)
    {
        fprintf(stderr, "Embree init failed: rtcNewDevice returned null.\n");
        std::abort();
    }
}

CPUDevice::~CPUDevice()
{
    if (embreeDevice)
    {
        rtcReleaseDevice(embreeDevice);
        embreeDevice = nullptr;
    }
}

DeviceKind CPUDevice::GetKind() const
{
    return DeviceKind::CPU;
}

DeviceMemoryView<uint8_t> CPUDevice::AllocBytes(size_t numBytes)
{
    YBI_ASSERT(numBytes != 0);
    uint8_t *ptr = (uint8_t *)std::malloc(numBytes);
    YBI_ASSERT(ptr);
    totalAllocated += numBytes;
    return {ptr, numBytes};
}

void CPUDevice::FreeBytes(DeviceMemoryView<uint8_t> &view)
{
    if (view.data() == nullptr || view.size() == 0)
    {
        return;
    }
    totalAllocated -= view.numBytes();
    std::free(view.data());
    view = {};
}

void CPUDevice::CopyBytesToDevice(DeviceMemoryView<uint8_t> dst,
                                  const void *src,
                                  size_t numBytes)
{
    YBI_ASSERT(dst.data());
    YBI_ASSERT(src);
    YBI_ASSERT(numBytes <= dst.numBytes());
    std::memcpy(dst.data(), src, numBytes);
}

void CPUDevice::CopyBytesToHost(void *dst,
                                DeviceMemoryView<const uint8_t> src,
                                size_t numBytes)
{
    YBI_ASSERT(dst);
    YBI_ASSERT(src.data());
    YBI_ASSERT(numBytes <= src.numBytes());
    std::memcpy(dst, src.data(), numBytes);
}

size_t CPUDevice::GetBVHAllocatedBytes() const
{
    return bvhTotalAllocated;
}

void CPUDevice::BuildBVH(Scene *scene)
{
    BuildEmbreeBVH(this, scene);
}

YBI_NAMESPACE_END

#endif

