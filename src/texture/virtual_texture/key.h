#pragma once

#include <math.h>
#include <stdint.h>

#include "device/util.h"
#include "util/math_common.h"

namespace ybi
{
namespace texture
{

YBI_DEVICE unsigned long long PackVirtualTextureKey(unsigned int tileX,
                                                    unsigned int tileY,
                                                    unsigned int udimBits,
                                                    unsigned int textureId,
                                                    unsigned int mip)
{
    return (static_cast<unsigned long long>(tileX & 0x1ffu) << 0u) |
           (static_cast<unsigned long long>(tileY & 0x1ffu) << 9u) |
           (static_cast<unsigned long long>(udimBits & 0x7fu) << 18u) |
           (static_cast<unsigned long long>(textureId & 0x7fffffu) << 25u) |
           (static_cast<unsigned long long>(mip & 0xfu) << 48u);
}

YBI_DEVICE unsigned int UnpackVirtualTextureKeyTileX(unsigned long long key)
{
    return static_cast<unsigned int>((key >> 0u) & 0x1ffull);
}

YBI_DEVICE unsigned int UnpackVirtualTextureKeyTileY(unsigned long long key)
{
    return static_cast<unsigned int>((key >> 9u) & 0x1ffull);
}

YBI_DEVICE unsigned int UnpackVirtualTextureKeyUdim(unsigned long long key)
{
    return 1001u + static_cast<unsigned int>((key >> 18u) & 0x7full);
}

YBI_DEVICE unsigned int UnpackVirtualTextureKeyTextureId(unsigned long long key)
{
    return static_cast<unsigned int>((key >> 25u) & 0x7fffffull);
}

YBI_DEVICE unsigned int UnpackVirtualTextureKeyMip(unsigned long long key)
{
    return static_cast<unsigned int>((key >> 48u) & 0xfull);
}

YBI_DEVICE int UdimFromUV(float u, float v)
{
    const int udimU = int(floorf(u));
    const int udimV = int(floorf(v));
    return 1001 + udimU + 10 * udimV;
}

YBI_DEVICE unsigned int UdimBitsFromUdim(int udim)
{
    return static_cast<unsigned int>(Clamp(udim - 1001, 0, 127));
}

YBI_DEVICE int UdimSlotFromUdim(int udim, int stride)
{
    return Clamp(udim - 1001, 0, Max(stride - 1, 0));
}

YBI_DEVICE float UnitUV(float v)
{
    return v - floorf(v);
}

YBI_DEVICE int TexelFromUnitUV(float unit, int textureSize)
{
    const int maxTexel = Max(textureSize - 1, 0);
    return Clamp(int(floorf(unit * float(textureSize))), 0, maxTexel);
}

YBI_DEVICE unsigned int TileCoordFromTexel(int texel, int tileSize)
{
    const int safeTileSize = Max(tileSize, 1);
    return static_cast<unsigned int>(Clamp(texel / safeTileSize, 0, 511));
}

} // namespace texture
} // namespace ybi
