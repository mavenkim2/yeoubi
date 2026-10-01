#pragma once

#include "core/assert.h"
#include "core/math.h"
#include "core/vec2f.h"
#include "core/vec4f.h"
#include "device/cuda_compat.h"
#include "device/device_texture.h"
#include "device/kernel_params.h"
#include "render/texture_sampling.h"
#include "render/virtual_texture_interop.h"
#include "render/virtual_texture_sampling.h"
#include <cuda_runtime.h>

namespace ybi
{

YBI_GPU uint32_t CalculateMaxVirtualTextureMip(DeviceTextureDescriptor descriptor,
                                               uint32_t pageShift)
{
    const uint32_t pageSize = 1u << pageShift;
    uint32_t maxMip = 0;

    while (maxMip < 15u)
    {
        const uint32_t nextMip = maxMip + 1u;
        if ((descriptor.width >> nextMip) < pageSize || (descriptor.height >> nextMip) < pageSize)
        {
            break;
        }

        maxMip = nextMip;
    }

    return maxMip;
}

YBI_GPU Vec4f SamplePointVirtualTexture(const KernelParams &params,
                                        DeviceTextureDescriptor descriptor,
                                        uint32_t pageShift,
                                        Vec2f uv,
                                        uint32_t mip,
                                        UnpackedFeedbackRequest *request = nullptr,
                                        uint32_t udimTableIndex = 0,
                                        uint32_t base = 0,
                                        uint32_t offset = 0)
{
    const uint32_t pageSize = (1u << pageShift);
    const uint32_t paddedPageSize = pageSize + 2 * descriptor.gutterSize;
    Vec2f normalizedUv = uv - Floor(uv);
    // Floating-point cancellation can make fract(uv) exactly 1 for values just below an
    // integer, particularly for negative or large UVs. Keep the repeated coordinate in
    // [0, 1) so it cannot produce a one-past-the-end texel or page.
    constexpr float maxNormalizedUv = 0x1.fffffep-1f;
    normalizedUv.x = Clamp(normalizedUv.x, 0.f, maxNormalizedUv);
    normalizedUv.y = Clamp(normalizedUv.y, 0.f, maxNormalizedUv);

    uint32_t mipWidth = descriptor.width >> mip;
    uint32_t mipHeight = descriptor.height >> mip;

    if (request)
    {
        request->textureIndex = descriptor.textureIndex;
        request->tileX = 0;
        request->tileY = 0;
        request->mip = mip;
        request->pageDescriptorIndex = ~0u;
    }

    if (mipWidth < pageSize || mipHeight < pageSize)
    {
        return descriptor.fallback;
    }

    Vec2f texel = normalizedUv * Vec2f(mipWidth, mipHeight);
    const Vec2f page = Floor(texel / float(pageSize));
    const Vec2f texelInPage = texel - page * pageSize + descriptor.gutterSize;

    const Vec2f pageUv = texelInPage / paddedPageSize;

    const uint32_t pageOffset = CalculatePageDescriptorOffset(descriptor.width,
                                                              descriptor.height,
                                                              mip,
                                                              pageShift,
                                                              static_cast<uint32_t>(page.x),
                                                              static_cast<uint32_t>(page.y));
    if (pageOffset >= descriptor.pageDescriptorCount)
    {
        printf("%u %u, %u: %u %u %f %f %f %f, %u %u\n",
               pageOffset,
               descriptor.pageDescriptorCount,
               udimTableIndex,
               descriptor.width,
               descriptor.height,
               page.x,
               page.y,
               uv.x,
               uv.y,
               base,
               offset);
        YBI_ASSERT(pageOffset < descriptor.pageDescriptorCount);
    }
    const uint32_t pageIndex = descriptor.pageDescriptorOffset + pageOffset;

    if (request)
    {
        request->tileX = static_cast<uint32_t>(page.x);
        request->tileY = static_cast<uint32_t>(page.y);
        request->pageDescriptorIndex = pageIndex;
    }

    const PageDescriptor pageDescriptor = params.pageDescriptors[pageIndex];

    if (pageDescriptor.bits == ~0u)
    {
        return descriptor.fallback;
    }

    const uint32_t layer = pageDescriptor.bits >> 27u;
    const uint32_t textureIndex = pageDescriptor.bits & 0x07ffffffu;

    // DeviceTexture texture = params.texturePools[textureIndex];
    uint64_t textureObject = pageDescriptor.textureObject;
    if (textureObject == 0ull)
    {
        return descriptor.fallback;
    }

    Vec4f result =
        SampleDeviceTexture(textureObject, descriptor.format, pageUv.x, pageUv.y, layer);

    return result;
}

YBI_GPU void WriteUniqueFeedback(const KernelParams &params, UnpackedFeedbackRequest request)
{
    const uint32_t pageDescriptorIndex = request.pageDescriptorIndex;
    if (pageDescriptorIndex == ~0u)
    {
        return;
    }

    uint32_t *bits = params.pageDescriptorBits;

    const uint32_t wordIndex = pageDescriptorIndex >> 5u;
    const uint32_t bitIndex = pageDescriptorIndex & 31u;

    uint32_t old = atomicOr(&bits[wordIndex], 1u << bitIndex);
    if (((old >> bitIndex) & 1u) == 0)
    {
        uint32_t requestIndex = atomicAdd(&params.feedbackRequests[0].textureIndex, 1u);

        params.feedbackRequests[requestIndex + 1] =
            PackFeedbackRequest(request.textureIndex, request.tileX, request.tileY, request.mip);
    }
}

} // namespace ybi
