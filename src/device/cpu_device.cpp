#include "device/cpu_device.h"
#include "render/dispatch_types.h"

#if defined(WITH_EMBREE)

#include "util/aligned_malloc.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>

YBI_NAMESPACE_BEGIN

namespace
{
static size_t KernelIndex(RenderKernelId kernel)
{
    return static_cast<size_t>(kernel);
}

static bool CPUPrimaryDiffuseStub(CPUDevice *device, const DispatchParams &params)
{
    (void)device;
    (void)params;
    fprintf(stderr, "CPU primary diffuse kernel is not wired yet.\n");
    return false;
}

static bool CPUAOStub(CPUDevice *device, const DispatchParams &params)
{
    (void)device;
    (void)params;
    fprintf(stderr, "CPU AO kernel is not wired yet.\n");
    return false;
}
} // namespace

CPUDevice::CPUDevice()
{
    embreeDevice = rtcNewDevice(nullptr);
    if (!embreeDevice)
    {
        fprintf(stderr, "Embree init failed: rtcNewDevice returned null.\n");
        std::abort();
    }
    RegisterKernel(RenderKernelId::PrimaryDiffuse, CPUPrimaryDiffuseStub);
    RegisterKernel(RenderKernelId::AO, CPUAOStub);
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
    constexpr size_t kCPUDeviceAlignment = 64;
    uint8_t *ptr = (uint8_t *)util::AlignedAlloc(numBytes, kCPUDeviceAlignment);
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
    util::AlignedFree(view.data());
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

bool CPUDevice::CreateTexture(const DeviceTextureCreateInfo &info,
                              DeviceTexture *outTexture,
                              std::string *outError)
{
    if (!outTexture)
    {
        if (outError)
        {
            *outError = "CreateTexture: outTexture is null";
        }
        return false;
    }
    *outTexture = {};

    if (info.format != DeviceTextureFormat::RGBA8_UNORM)
    {
        if (outError)
        {
            *outError = "CreateTexture: unsupported format";
        }
        return false;
    }
    if (info.pixels == nullptr || info.width == 0 || info.height == 0)
    {
        if (outError)
        {
            *outError = "CreateTexture: invalid texture input";
        }
        return false;
    }

    const size_t expectedBytes = size_t(info.width) * size_t(info.height) * 4u;
    if (info.pixelBytes < expectedBytes)
    {
        if (outError)
        {
            *outError = "CreateTexture: pixelBytes too small";
        }
        return false;
    }

    constexpr size_t kTextureAlignment = 64;
    void *storage = util::AlignedAlloc(expectedBytes, kTextureAlignment);
    if (!storage)
    {
        if (outError)
        {
            *outError = "CreateTexture: allocation failed";
        }
        return false;
    }

    std::memcpy(storage, info.pixels, expectedBytes);

    outTexture->handle = reinterpret_cast<uint64_t>(storage);
    outTexture->allocation = reinterpret_cast<uint64_t>(storage);
    outTexture->width = info.width;
    outTexture->height = info.height;
    outTexture->wrapS = info.wrapS;
    outTexture->wrapT = info.wrapT;
    outTexture->filter = info.filter;
    outTexture->format = info.format;
    outTexture->valid = true;
    return true;
}

void CPUDevice::DestroyTexture(DeviceTexture &texture)
{
    if (texture.allocation != 0)
    {
        util::AlignedFree(reinterpret_cast<void *>(texture.allocation));
    }
    texture = {};
}

size_t CPUDevice::GetBVHAllocatedBytes() const
{
    return bvhTotalAllocated;
}

void CPUDevice::BuildBVH(Scene *scene)
{
    BuildEmbreeBVH(this, scene);
}

bool CPUDevice::DispatchKernel(RenderKernelId kernel, const DispatchParams &params)
{
    const size_t idx = KernelIndex(kernel);
    if (idx >= kernels.size() || kernels[idx] == nullptr)
    {
        fprintf(stderr, "CPUDevice::DispatchKernel missing kernel %zu.\n", idx);
        return false;
    }
    return kernels[idx](this, params);
}

void CPUDevice::RegisterKernel(RenderKernelId kernel, KernelFn fn)
{
    const size_t idx = KernelIndex(kernel);
    YBI_ASSERT(idx < kernels.size());
    kernels[idx] = fn;
}

YBI_NAMESPACE_END

#endif
