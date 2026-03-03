#pragma once

#if defined(WITH_EMBREE)

#include "device/device.h"
#include "util/assert.h"
#include "util/base.h"
#include <cstddef>

#if __has_include(<embree4/rtcore.h>)
#include <embree4/rtcore.h>
#elif __has_include(<embree3/rtcore.h>)
#include <embree3/rtcore.h>
#else
#error "Embree headers not found: expected embree3/rtcore.h or embree4/rtcore.h"
#endif

YBI_NAMESPACE_BEGIN

struct Scene;

struct CPUDevice : Device
{
    RTCDevice embreeDevice = nullptr;
    size_t totalAllocated = 0;
    size_t bvhTotalAllocated = 0;

    CPUDevice();
    ~CPUDevice() override;

    DeviceKind GetKind() const override;
    DeviceMemoryView<uint8_t> AllocBytes(size_t numBytes) override;
    void FreeBytes(DeviceMemoryView<uint8_t> &view) override;
    void CopyBytesToDevice(DeviceMemoryView<uint8_t> dst,
                           const void *src,
                           size_t numBytes) override;
    void CopyBytesToHost(void *dst,
                         DeviceMemoryView<const uint8_t> src,
                         size_t numBytes) override;
    size_t GetBVHAllocatedBytes() const override;
    void BuildBVH(Scene *scene) override;
    bool DispatchKernel(RenderKernelId kernel, const DispatchParams &params) override;
};

void BuildEmbreeBVH(CPUDevice *cpuDevice, Scene *scene);

YBI_NAMESPACE_END

#endif
