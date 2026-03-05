#pragma once

#include <stdint.h>
#include <math.h>

#if defined(__CUDACC__)
#    define YBI_VT_HD __host__ __device__ __forceinline__
#else
#    define YBI_VT_HD inline
#endif

namespace ybi
{
namespace texture
{

YBI_VT_HD int VtClampInt(int v, int lo, int hi)
{
    return v < lo ? lo : (v > hi ? hi : v);
}

YBI_VT_HD int VtMaxInt(int a, int b)
{
    return a > b ? a : b;
}

YBI_VT_HD unsigned long long PackVirtualTextureKey(unsigned int tileX,
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

YBI_VT_HD unsigned int UnpackVirtualTextureKeyTileX(unsigned long long key)
{
    return static_cast<unsigned int>((key >> 0u) & 0x1ffull);
}

YBI_VT_HD unsigned int UnpackVirtualTextureKeyTileY(unsigned long long key)
{
    return static_cast<unsigned int>((key >> 9u) & 0x1ffull);
}

YBI_VT_HD unsigned int UnpackVirtualTextureKeyUdim(unsigned long long key)
{
    return 1001u + static_cast<unsigned int>((key >> 18u) & 0x7full);
}

YBI_VT_HD unsigned int UnpackVirtualTextureKeyTextureId(unsigned long long key)
{
    return static_cast<unsigned int>((key >> 25u) & 0x7fffffull);
}

YBI_VT_HD unsigned int UnpackVirtualTextureKeyMip(unsigned long long key)
{
    return static_cast<unsigned int>((key >> 48u) & 0xfull);
}

YBI_VT_HD int UdimFromUV(float u, float v)
{
    const int udimU = int(floorf(u));
    const int udimV = int(floorf(v));
    return 1001 + udimU + 10 * udimV;
}

YBI_VT_HD unsigned int UdimBitsFromUdim(int udim)
{
    return static_cast<unsigned int>(VtClampInt(udim - 1001, 0, 127));
}

YBI_VT_HD int UdimSlotFromUdim(int udim, int stride)
{
    return VtClampInt(udim - 1001, 0, VtMaxInt(stride - 1, 0));
}

YBI_VT_HD float UnitUV(float v)
{
    return v - floorf(v);
}

YBI_VT_HD int TexelFromUnitUV(float unit, int textureSize)
{
    const int maxTexel = VtMaxInt(textureSize - 1, 0);
    return VtClampInt(int(floorf(unit * float(textureSize))), 0, maxTexel);
}

YBI_VT_HD unsigned int TileCoordFromTexel(int texel, int tileSize)
{
    const int safeTileSize = VtMaxInt(tileSize, 1);
    return static_cast<unsigned int>(VtClampInt(texel / safeTileSize, 0, 511));
}

} // namespace texture
} // namespace ybi

#undef YBI_VT_HD
