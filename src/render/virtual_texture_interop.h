#pragma once

#include "core/assert.h"
#include "core/bit.h"
#include "core/vec4f.h"
#include "device/cuda_compat.h"
#include "device/device_texture.h"
#include <cassert>
#include <cstdint>

namespace ybi
{

struct MaterialTextureHandle
{
    uint32_t udimBase = 0;
    uint32_t udimTableOffset = 0;
};

struct PageDescriptorUpdateRequest
{
    uint32_t pageDescriptorIndex;
    uint64_t textureObject;
    uint32_t bits;
};

struct PageDescriptor
{
    uint64_t textureObject;
    uint32_t bits;
};

struct DeviceTextureDescriptor
{
    uint32_t textureIndex;
    TextureFormat format;
    uint32_t pageDescriptorOffset;
    uint32_t pageDescriptorCount;
    uint32_t width;
    uint32_t height;
    Vec4f fallback;
};

struct UnpackedFeedbackRequest
{
    uint32_t textureIndex;
    uint32_t tileX;
    uint32_t tileY;
    uint32_t mip;
    uint32_t pageDescriptorIndex;
};

struct FeedbackRequest
{
    uint32_t textureIndex;
    // 8, 8, 4 bits
    uint32_t tileX_tileY_mip;
};

YBI_CPU_GPU uint32_t CalculateVirtualTexturePageCount(uint32_t dimension, uint32_t pageShift)
{
    YBI_ASSERT(pageShift < 32u);
    YBI_ASSERT(dimension > 0u);

    if (dimension == 0u)
    {
        return 0u;
    }

    // Ceiling division by the power-of-two page size without overflowing dimension.
    return (dimension + (1u << pageShift) - 1u) >> pageShift;
}

YBI_CPU_GPU uint32_t CalculatePageDescriptorOffset(uint32_t width,
                                                   uint32_t height,
                                                   uint32_t mip,
                                                   uint32_t pageShift,
                                                   uint32_t pageX,
                                                   uint32_t pageY)
{
    YBI_ASSERT(width > 0u && height > 0u);

    uint32_t offset = 0;
    for (uint32_t i = 0; i < mip; i++)
    {
        const uint32_t pagesX = CalculateVirtualTexturePageCount(width, pageShift);
        const uint32_t pagesY = CalculateVirtualTexturePageCount(height, pageShift);
        offset += pagesX * pagesY;
        width = width > 1u ? width >> 1u : 1u;
        height = height > 1u ? height >> 1u : 1u;
    }

    const uint32_t pagesX = CalculateVirtualTexturePageCount(width, pageShift);
    const uint32_t pagesY = CalculateVirtualTexturePageCount(height, pageShift);
    YBI_ASSERT(pageX < pagesX);
    YBI_ASSERT(pageY < pagesY);

    return offset + pageY * pagesX + pageX;
}

YBI_CPU_GPU FeedbackRequest PackFeedbackRequest(uint32_t textureIndex,
                                                uint32_t tileX,
                                                uint32_t tileY,
                                                uint32_t mip)
{
    assert(tileX < 256);
    assert(tileY < 256);
    assert(mip < 16);

    FeedbackRequest request;
    request.textureIndex = textureIndex;
    request.tileX_tileY_mip = tileX | (tileY << 8u) | (mip << 16u);
    return request;
}

YBI_CPU_GPU void UnpackFeedbackRequest(FeedbackRequest request,
                                       uint32_t &textureIndex,
                                       uint32_t &tileX,
                                       uint32_t &tileY,
                                       uint32_t &mip)
{
    textureIndex = request.textureIndex;
    tileX = request.tileX_tileY_mip & 0xffu;
    tileY = (request.tileX_tileY_mip >> 8u) & 0xffu;
    mip = (request.tileX_tileY_mip >> 16u) & 0xfu;
}

} // namespace ybi
