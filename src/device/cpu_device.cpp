#include "device/cpu_device.h"
#include "render/cpu_kernels.h"
#if defined(WITH_EMBREE)

#include "util/aligned_malloc.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>

namespace ybi
{

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

bool CPUDevice::InitializeKernels(const std::string &, std::string *)
{
    return true;
}

void CPUDevice::DestroyKernels()
{
}

void CPUDevice::ClearTransientMemory()
{
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

    const size_t pixelBytes = TextureFormatPixelBytes(info.format);
    if (pixelBytes == 0u)
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

    const size_t expectedBytes = size_t(info.width) * size_t(info.height) * pixelBytes;
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

bool CPUDevice::UpdateTextureRegion(const DeviceTexture &texture,
                                    uint32_t x,
                                    uint32_t y,
                                    uint32_t width,
                                    uint32_t height,
                                    const void *pixels,
                                    size_t pixelBytes,
                                    std::string *outError)
{
    if (!texture.valid || texture.handle == 0 || pixels == nullptr)
    {
        if (outError)
        {
            *outError = "UpdateTextureRegion: invalid texture input";
        }
        return false;
    }

    const size_t texelBytes = TextureFormatPixelBytes(texture.format);
    if (texelBytes == 0u || width == 0u || height == 0u || x + width > texture.width ||
        y + height > texture.height)
    {
        if (outError)
        {
            *outError = "UpdateTextureRegion: invalid region";
        }
        return false;
    }

    const size_t srcRowBytes = static_cast<size_t>(width) * texelBytes;
    const size_t expectedBytes = srcRowBytes * static_cast<size_t>(height);
    if (pixelBytes < expectedBytes)
    {
        if (outError)
        {
            *outError = "UpdateTextureRegion: pixelBytes too small";
        }
        return false;
    }

    uint8_t *dstPixels = reinterpret_cast<uint8_t *>(texture.handle);
    const size_t dstRowBytes = static_cast<size_t>(texture.width) * texelBytes;
    const size_t dstBase =
        (static_cast<size_t>(y) * static_cast<size_t>(texture.width) + static_cast<size_t>(x)) *
        texelBytes;
    const uint8_t *srcPixels = reinterpret_cast<const uint8_t *>(pixels);
    for (uint32_t row = 0u; row < height; ++row)
    {
        std::memcpy(dstPixels + dstBase + static_cast<size_t>(row) * dstRowBytes,
                    srcPixels + static_cast<size_t>(row) * srcRowBytes,
                    srcRowBytes);
    }
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
    return CPUDispatchKernel(params, kernel);
}

} // namespace ybi

#endif
