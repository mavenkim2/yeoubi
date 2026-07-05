#pragma once

#include "core/assert.h"
#include "core/bit.h"
#include "core/vec4f.h"
#include "device/cuda_compat.h"
#include <cstdint>

namespace ybi
{

struct PageDescriptorUpdateRequest
{
    uint32_t pageDescriptorIndex;
    uint32_t bits;
};

struct PageDescriptor
{
    // texture, page
    uint32_t bits;
};

struct DeviceTextureDescriptor
{
    uint32_t pageDescriptorOffset;
    uint32_t width;
    uint32_t height;
    Vec4f fallback;
};

YBI_CPU_GPU uint32_t CalculatePageDescriptorOffset(uint32_t width,
                                                   uint32_t height,
                                                   uint32_t mip,
                                                   uint32_t pageShift,
                                                   uint32_t pageX,
                                                   uint32_t pageY)
{
    YBI_ASSERT(IsPowerOfTwo(width) && IsPowerOfTwo(height));

    uint32_t basePagesX = width >> pageShift;
    uint32_t basePagesY = height >> pageShift;
    uint32_t basePages = basePagesX * basePagesY;

    uint32_t mipPagesX = basePagesX >> mip;
    uint32_t mipPagesBefore = ((basePages - (basePages >> (2 * mip))) * 4) / 3;

    uint32_t offset = mipPagesBefore + pageY * mipPagesX + pageX;
    return offset;
}

} // namespace ybi
