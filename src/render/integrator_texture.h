#pragma once

#include "render/integrator_common.h"
#include "render/launch_params.h"
#include "texture/virtual_texture/key.h"

namespace ybi
{
namespace render
{
namespace integrator
{

YBI_INTEGRATOR_HD Vec3 MaterialSampleToViewColor(const LaunchParams &params, const Vec4 &sample)
{
    if (params.textureViewSemantic == kSemanticRoughness ||
        params.textureViewSemantic == kSemanticMetallic ||
        params.textureViewSemantic == kSemanticOcclusion ||
        params.textureViewSemantic == kSemanticIor ||
        params.textureViewSemantic == kSemanticOpacity)
    {
        return MakeVec3(sample.x, sample.x, sample.x);
    }
    return MakeVec3(sample.x, sample.y, sample.z);
}

YBI_INTEGRATOR_HD uint32_t PackVirtualTexturePageEntry(unsigned int pageX,
                                                       unsigned int pageY,
                                                       unsigned int pageType,
                                                       unsigned int flags)
{
    return ((pageX & 0xfffu) << 0u) | ((pageY & 0xfffu) << 12u) | ((pageType & 0xfu) << 24u) |
           ((flags & 0xfu) << 28u);
}

YBI_INTEGRATOR_HD void UnpackVirtualTexturePageEntry(uint32_t packed,
                                                     unsigned int *outPageX,
                                                     unsigned int *outPageY,
                                                     unsigned int *outPageType,
                                                     unsigned int *outFlags)
{
    if (outPageX)
    {
        *outPageX = (packed >> 0u) & 0xfffu;
    }
    if (outPageY)
    {
        *outPageY = (packed >> 12u) & 0xfffu;
    }
    if (outPageType)
    {
        *outPageType = (packed >> 24u) & 0xfu;
    }
    if (outFlags)
    {
        *outFlags = (packed >> 28u) & 0xfu;
    }
}

YBI_INTEGRATOR_HD bool ResolveVirtualTextureInfo(const LaunchParams &params,
                                                 const LaunchParams::InstanceGeomRef &geomRef,
                                                 unsigned int mip,
                                                 unsigned int udimBits,
                                                 unsigned int tileX,
                                                 unsigned int tileY,
                                                 const LaunchParams::VirtualTextureTextureMeta **outMeta,
                                                 const LaunchParams::VirtualTextureMipInfo **outMipInfo,
                                                 unsigned int *outLocalUdim,
                                                 unsigned int *outVax,
                                                 unsigned int *outVay)
{
    if (params.virtualTextureTextureMeta == 0ull || params.virtualTextureMipInfos == 0ull ||
        geomRef.materialIndex < 0 || geomRef.materialIndex >= params.virtualTextureTextureMetaCount)
    {
        return false;
    }

    const LaunchParams::VirtualTextureTextureMeta *textures =
        reinterpret_cast<const LaunchParams::VirtualTextureTextureMeta *>(params.virtualTextureTextureMeta);
    const LaunchParams::VirtualTextureTextureMeta &meta = textures[geomRef.materialIndex];
    if (mip >= meta.mipCount || udimBits >= 128u)
    {
        return false;
    }
    const int localUdim = static_cast<int>(meta.udimToLocal[udimBits]);
    if (localUdim < 0)
    {
        return false;
    }

    const LaunchParams::VirtualTextureMipInfo *mipInfos =
        reinterpret_cast<const LaunchParams::VirtualTextureMipInfo *>(params.virtualTextureMipInfos);
    const LaunchParams::VirtualTextureMipInfo &mipInfo = mipInfos[meta.mipInfoOffset + mip];
    if (tileX >= mipInfo.pagesPerUdimX || tileY >= mipInfo.pagesPerUdimY)
    {
        return false;
    }

    const unsigned int vaX =
        mipInfo.basePageX + static_cast<unsigned int>(localUdim) * mipInfo.pagesPerUdimX + tileX;
    const unsigned int vaY = mipInfo.basePageY + tileY;
    if (vaX >= mipInfo.basePageX + mipInfo.pageCountX || vaY >= mipInfo.basePageY + mipInfo.pageCountY)
    {
        return false;
    }

    if (outMeta)
    {
        *outMeta = &meta;
    }
    if (outMipInfo)
    {
        *outMipInfo = &mipInfo;
    }
    if (outLocalUdim)
    {
        *outLocalUdim = static_cast<unsigned int>(localUdim);
    }
    if (outVax)
    {
        *outVax = vaX;
    }
    if (outVay)
    {
        *outVay = vaY;
    }
    return true;
}

YBI_INTEGRATOR_HD bool ReadVirtualTexturePageEntry(const LaunchParams &params,
                                                   unsigned int mip,
                                                   unsigned int vaX,
                                                   unsigned int vaY,
                                                   uint32_t *outEntry)
{
    if (!outEntry || params.virtualTexturePageTableEntries == 0ull ||
        params.virtualTexturePageTableMipOffsets == 0ull ||
        params.virtualTexturePageTableMipWidths == 0ull ||
        params.virtualTexturePageTableMipHeights == 0ull || mip >= uint32_t(params.virtualTexturePageTableMipCount))
    {
        return false;
    }

    const uint32_t *offsets =
        reinterpret_cast<const uint32_t *>(params.virtualTexturePageTableMipOffsets);
    const uint32_t *widths =
        reinterpret_cast<const uint32_t *>(params.virtualTexturePageTableMipWidths);
    const uint32_t *heights =
        reinterpret_cast<const uint32_t *>(params.virtualTexturePageTableMipHeights);
    const uint32_t width = widths[mip];
    const uint32_t height = heights[mip];
    if (vaX >= width || vaY >= height)
    {
        return false;
    }

    const uint32_t *entries =
        reinterpret_cast<const uint32_t *>(params.virtualTexturePageTableEntries);
    *outEntry = entries[offsets[mip] + vaY * width + vaX];
    return true;
}

template <typename State>
YBI_INTEGRATOR_HD void TryWriteTextureFeedback(State &state,
                                               const LaunchParams::InstanceGeomRef &geomRef,
                                               unsigned int primitiveIndex,
                                               float u,
                                               float v,
                                               float uu,
                                               float vv,
                                               int textureWidth,
                                               int textureHeight,
                                               unsigned int mip)
{
    const LaunchParams &params = state.Params();
    if (params.feedbackKeys == 0ull || params.feedbackStats == 0ull || params.feedbackCapacity <= 0 ||
        params.feedbackSamplePercent <= 0)
    {
        return;
    }

    const UInt2 launchIndex = state.LaunchIndex();
    unsigned int seed =
        Hash32((launchIndex.x + 1u) * 73856093u ^ (launchIndex.y + 1u) * 19349663u ^
               (primitiveIndex + 1u) * 83492791u ^
               (static_cast<unsigned int>(params.currentSpp) + 1u) * 2654435761u ^
               (static_cast<unsigned int>(ClampInt(geomRef.materialIndex, 0, int((1u << 23u) - 1u))) +
                1u) *
                   374761393u);
    if ((seed % 100u) >= static_cast<unsigned int>(params.feedbackSamplePercent))
    {
        return;
    }

    const int tileSize = MaxInt(params.feedbackTileSize, 1);
    const int mipInt = ClampInt(int(mip), 0, 30);
    const int mipWidth = MaxInt(textureWidth >> mipInt, 1);
    const int mipHeight = MaxInt(textureHeight >> mipInt, 1);
    const int texelX = ybi::texture::TexelFromUnitUV(uu, mipWidth);
    const int texelY = ybi::texture::TexelFromUnitUV(vv, mipHeight);
    const unsigned int tileX = ybi::texture::TileCoordFromTexel(texelX, tileSize);
    const unsigned int tileY = ybi::texture::TileCoordFromTexel(texelY, tileSize);

    const int udim = ybi::texture::UdimFromUV(u, v);
    const unsigned int udimBits = ybi::texture::UdimBitsFromUdim(udim);

    const unsigned int textureId =
        static_cast<unsigned int>(ClampInt(geomRef.materialIndex, 0, int((1u << 23u) - 1u)));
    const unsigned long long key =
        ybi::texture::PackVirtualTextureKey(tileX, tileY, udimBits, textureId, mip);

    state.RecordFeedbackKey(key);
}

template <typename State>
YBI_INTEGRATOR_HD bool TrySampleVirtualTexture(State &state,
                                               const LaunchParams::InstanceGeomRef &geomRef,
                                               unsigned int primitiveIndex,
                                               const LaunchParams::MaterialTextureRef &textureRef,
                                               float u,
                                               float v,
                                               float uu,
                                               float vv,
                                               Vec3 &outColor)
{
    const LaunchParams &params = state.Params();
    if (params.virtualTextureEnabled == 0)
    {
        return false;
    }

    if (params.virtualTexturePageTableEntries == 0ull || params.virtualTextureStreamPixels == 0ull ||
        params.virtualTextureTextureMeta == 0ull || params.virtualTextureMipInfos == 0ull)
    {
        outColor = MakeVec3(1.0f, 0.0f, 1.0f);
        return true;
    }

    bool blackS = false;
    bool blackT = false;
    const float wrappedU = ApplyWrapMode(uu, textureRef.wrapS, blackS);
    const float wrappedV = ApplyWrapMode(vv, textureRef.wrapT, blackT);
    if (blackS || blackT)
    {
        outColor = MakeVec3(0.0f, 0.0f, 0.0f);
        return true;
    }

    const unsigned int mip = static_cast<unsigned int>(MaxInt(params.virtualTextureSampleMip, 0));
    const int mipWidth = MaxInt(textureRef.width >> int(mip), 1);
    const int mipHeight = MaxInt(textureRef.height >> int(mip), 1);
    const int tileSize = MaxInt(params.virtualTexturePageSize, 1);
    const int texelX = ybi::texture::TexelFromUnitUV(wrappedU, mipWidth);
    const int texelY = ybi::texture::TexelFromUnitUV(wrappedV, mipHeight);
    const unsigned int tileX = ybi::texture::TileCoordFromTexel(texelX, tileSize);
    const unsigned int tileY = ybi::texture::TileCoordFromTexel(texelY, tileSize);
    const int udim = ybi::texture::UdimFromUV(u, v);
    const unsigned int udimBits = ybi::texture::UdimBitsFromUdim(udim);
    const unsigned int textureId =
        static_cast<unsigned int>(ClampInt(geomRef.materialIndex, 0, int((1u << 23u) - 1u)));
    const unsigned long long key =
        ybi::texture::PackVirtualTextureKey(tileX, tileY, udimBits, textureId, mip);

    const LaunchParams::VirtualTextureTextureMeta *meta = nullptr;
    const LaunchParams::VirtualTextureMipInfo *mipInfo = nullptr;
    unsigned int localUdim = 0u;
    unsigned int vaX = 0u;
    unsigned int vaY = 0u;
    if (!ResolveVirtualTextureInfo(
            params, geomRef, mip, udimBits, tileX, tileY, &meta, &mipInfo, &localUdim, &vaX, &vaY))
    {
        outColor = MakeVec3(1.0f, 0.0f, 1.0f);
        return true;
    }

    uint32_t packedEntry = 0u;
    if (!ReadVirtualTexturePageEntry(params, mipInfo->level, vaX, vaY, &packedEntry))
    {
        outColor = MakeVec3(1.0f, 0.0f, 1.0f);
        return true;
    }

    unsigned int pageX = 0u;
    unsigned int pageY = 0u;
    unsigned int pageType = 0u;
    unsigned int flags = 0u;
    UnpackVirtualTexturePageEntry(packedEntry, &pageX, &pageY, &pageType, &flags);

    const unsigned char *samplePixels = nullptr;
    unsigned long long sampleOffset = 0ull;
    const int localX = ClampInt(texelX - int(tileX) * tileSize, 0, tileSize - 1);
    const int localY = ClampInt(texelY - int(tileY) * tileSize, 0, tileSize - 1);
    if (pageType == kVirtualTexturePageTypeStream &&
        pageX < static_cast<unsigned int>(MaxInt(params.virtualTextureStreamPageCountX, 0)) &&
        pageY < static_cast<unsigned int>(MaxInt(params.virtualTextureStreamPageCountY, 0)))
    {
        const unsigned long long pageIndex =
            static_cast<unsigned long long>(pageY) *
                static_cast<unsigned long long>(params.virtualTextureStreamPageCountX) +
            static_cast<unsigned long long>(pageX);
        const unsigned long long pageBytes =
            static_cast<unsigned long long>(tileSize) * static_cast<unsigned long long>(tileSize) * 4ull;
        sampleOffset = pageIndex * pageBytes +
                       (static_cast<unsigned long long>(localY) * static_cast<unsigned long long>(tileSize) +
                        static_cast<unsigned long long>(localX)) *
                           4ull;
        samplePixels = reinterpret_cast<const unsigned char *>(params.virtualTextureStreamPixels);
    }
    else if (pageType == kVirtualTexturePageTypeTail && meta->tailPixels != 0ull &&
             pageX < meta->tailPageCountX && pageY < meta->tailPageCountY)
    {
        const int tailX = ybi::texture::TexelFromUnitUV(wrappedU, tileSize);
        const int tailY = ybi::texture::TexelFromUnitUV(wrappedV, tileSize);
        const unsigned long long pageIndex =
            static_cast<unsigned long long>(pageY) * static_cast<unsigned long long>(meta->tailPageCountX) +
            static_cast<unsigned long long>(pageX);
        const unsigned long long pageBytes =
            static_cast<unsigned long long>(tileSize) * static_cast<unsigned long long>(tileSize) * 4ull;
        sampleOffset = pageIndex * pageBytes +
                       (static_cast<unsigned long long>(tailY) * static_cast<unsigned long long>(tileSize) +
                        static_cast<unsigned long long>(tailX)) *
                           4ull;
        samplePixels = reinterpret_cast<const unsigned char *>(meta->tailPixels);
    }
    else if (meta->tailPixels != 0ull && localUdim < meta->activeUdimCount && meta->tailPageCountX > 0u)
    {
        const int tailX = ybi::texture::TexelFromUnitUV(wrappedU, tileSize);
        const int tailY = ybi::texture::TexelFromUnitUV(wrappedV, tileSize);
        const unsigned int fallbackPageX = localUdim % meta->tailPageCountX;
        const unsigned int fallbackPageY = localUdim / meta->tailPageCountX;
        const unsigned long long pageIndex =
            static_cast<unsigned long long>(fallbackPageY) * static_cast<unsigned long long>(meta->tailPageCountX) +
            static_cast<unsigned long long>(fallbackPageX);
        const unsigned long long pageBytes =
            static_cast<unsigned long long>(tileSize) * static_cast<unsigned long long>(tileSize) * 4ull;
        sampleOffset = pageIndex * pageBytes +
                       (static_cast<unsigned long long>(tailY) * static_cast<unsigned long long>(tileSize) +
                        static_cast<unsigned long long>(tailX)) *
                           4ull;
        samplePixels = reinterpret_cast<const unsigned char *>(meta->tailPixels);
    }
    else
    {
        outColor = MakeVec3(1.0f, 0.0f, 1.0f);
        return true;
    }

    const Vec4 sample = MakeVec4(float(samplePixels[sampleOffset + 0]) * (1.0f / 255.0f),
                                 float(samplePixels[sampleOffset + 1]) * (1.0f / 255.0f),
                                 float(samplePixels[sampleOffset + 2]) * (1.0f / 255.0f),
                                 float(samplePixels[sampleOffset + 3]) * (1.0f / 255.0f));
    outColor = MaterialSampleToViewColor(params, sample);
    TryWriteTextureFeedback(
        state, geomRef, primitiveIndex, u, v, wrappedU, wrappedV, textureRef.width, textureRef.height, mip);
    return true;
}

template <typename State>
YBI_INTEGRATOR_HD bool
TryWriteFeedbackOnly(State &state,
                     const LaunchParams::InstanceGeomRef &geomRef,
                     unsigned int primitiveIndex,
                     const Vec3 &barycentrics,
                     unsigned int mip)
{
    const LaunchParams &params = state.Params();
    if (params.materialTextureRefs == 0ull || params.materialTextureRefCount <= 0 ||
        params.materialTextureRefStride <= 0)
    {
        return false;
    }
    if (geomRef.texcoords == 0ull || geomRef.texcoordIndices == 0ull)
    {
        return false;
    }

    const int triCornerBase = int(primitiveIndex) * 3;
    if (triCornerBase + 2 >= geomRef.numTexcoordIndices)
    {
        return false;
    }

    const int *tcIndices = reinterpret_cast<const int *>(geomRef.texcoordIndices);
    const int t0 = tcIndices[triCornerBase + 0];
    const int t1 = tcIndices[triCornerBase + 1];
    const int t2 = tcIndices[triCornerBase + 2];
    if (t0 < 0 || t0 >= geomRef.numTexcoords || t1 < 0 || t1 >= geomRef.numTexcoords || t2 < 0 ||
        t2 >= geomRef.numTexcoords)
    {
        return false;
    }

    const UV2 *texcoords = reinterpret_cast<const UV2 *>(geomRef.texcoords);
    const UV2 uv0 = texcoords[t0];
    const UV2 uv1 = texcoords[t1];
    const UV2 uv2 = texcoords[t2];
    const float u = uv0.x * barycentrics.x + uv1.x * barycentrics.y + uv2.x * barycentrics.z;
    const float v = uv0.y * barycentrics.x + uv1.y * barycentrics.y + uv2.y * barycentrics.z;
    const int udim = ybi::texture::UdimFromUV(u, v);
    const int udimSlot = ybi::texture::UdimSlotFromUdim(udim, params.materialTextureRefStride);
    const float uu = ybi::texture::UnitUV(u);
    const float vv = ybi::texture::UnitUV(v);

    if (geomRef.materialIndex < 0 || geomRef.materialIndex >= params.materialTextureRefCount)
    {
        return false;
    }
    if (params.textureViewSemantic < 0 ||
        params.textureViewSemantic >= params.materialTextureRefSemanticCount)
    {
        return false;
    }

    const int materialIndex = geomRef.materialIndex;
    const LaunchParams::MaterialTextureRef *materialRefs =
        reinterpret_cast<const LaunchParams::MaterialTextureRef *>(params.materialTextureRefs);
    const int base = (materialIndex * params.materialTextureRefSemanticCount +
                      params.textureViewSemantic) *
                     params.materialTextureRefStride;
    const int slot = base + udimSlot;
    const int maxSlots = params.materialTextureRefCount * params.materialTextureRefSemanticCount *
                         params.materialTextureRefStride;
    if (slot < 0 || slot >= maxSlots)
    {
        return false;
    }

    const LaunchParams::MaterialTextureRef textureRef = materialRefs[slot];
    if (textureRef.width <= 0 || textureRef.height <= 0)
    {
        return false;
    }

    TryWriteTextureFeedback(
        state, geomRef, primitiveIndex, u, v, uu, vv, textureRef.width, textureRef.height, mip);
    return true;
}

template <typename State>
YBI_INTEGRATOR_HD bool TrySampleMaterialTexture(State &state,
                                                const LaunchParams::InstanceGeomRef &geomRef,
                                                unsigned int primitiveIndex,
                                                const Vec3 &barycentrics,
                                                Vec3 &outColor)
{
    const LaunchParams &params = state.Params();
    if (params.materialTextureRefs == 0ull || params.materialTextureRefCount <= 0 ||
        params.materialTextureRefStride <= 0)
    {
        return false;
    }
    if (geomRef.texcoords == 0ull || geomRef.texcoordIndices == 0ull)
    {
        return false;
    }

    const int triCornerBase = int(primitiveIndex) * 3;
    if (triCornerBase + 2 >= geomRef.numTexcoordIndices)
    {
        return false;
    }

    const int *tcIndices = reinterpret_cast<const int *>(geomRef.texcoordIndices);
    const int t0 = tcIndices[triCornerBase + 0];
    const int t1 = tcIndices[triCornerBase + 1];
    const int t2 = tcIndices[triCornerBase + 2];
    if (t0 < 0 || t0 >= geomRef.numTexcoords || t1 < 0 || t1 >= geomRef.numTexcoords || t2 < 0 ||
        t2 >= geomRef.numTexcoords)
    {
        return false;
    }

    const UV2 *texcoords = reinterpret_cast<const UV2 *>(geomRef.texcoords);
    const UV2 uv0 = texcoords[t0];
    const UV2 uv1 = texcoords[t1];
    const UV2 uv2 = texcoords[t2];
    const float u = uv0.x * barycentrics.x + uv1.x * barycentrics.y + uv2.x * barycentrics.z;
    const float v = uv0.y * barycentrics.x + uv1.y * barycentrics.y + uv2.y * barycentrics.z;
    const int udim = ybi::texture::UdimFromUV(u, v);
    const int udimSlot = ybi::texture::UdimSlotFromUdim(udim, params.materialTextureRefStride);
    const float uu = ybi::texture::UnitUV(u);
    const float vv = ybi::texture::UnitUV(v);

    if (geomRef.materialIndex < 0 || geomRef.materialIndex >= params.materialTextureRefCount)
    {
        return false;
    }
    const int materialIndex = geomRef.materialIndex;
    if (params.textureViewSemantic < 0 ||
        params.textureViewSemantic >= params.materialTextureRefSemanticCount)
    {
        return false;
    }

    const LaunchParams::MaterialTextureRef *materialRefs =
        reinterpret_cast<const LaunchParams::MaterialTextureRef *>(params.materialTextureRefs);
    const int base = (materialIndex * params.materialTextureRefSemanticCount +
                      params.textureViewSemantic) *
                     params.materialTextureRefStride;
    const int slot = base + udimSlot;
    const int maxSlots = params.materialTextureRefCount * params.materialTextureRefSemanticCount *
                         params.materialTextureRefStride;
    if (slot < 0 || slot >= maxSlots)
    {
        return false;
    }

    const LaunchParams::MaterialTextureRef textureRef = materialRefs[slot];
    if (params.virtualTextureEnabled != 0)
    {
        if (TrySampleVirtualTexture(state, geomRef, primitiveIndex, textureRef, u, v, uu, vv, outColor))
        {
            return true;
        }
    }

    if (textureRef.textureObject == 0ull || textureRef.valid == 0)
    {
        return false;
    }

    Vec4 sample = {};
    if (!state.SampleTexture2D(textureRef, uu, vv, sample))
    {
        return false;
    }
    outColor = MaterialSampleToViewColor(params, sample);
    TryWriteTextureFeedback(
        state, geomRef, primitiveIndex, u, v, uu, vv, textureRef.width, textureRef.height, 0u);
    state.MaybeLogSampleSuccess();
    return true;
}

} // namespace integrator
} // namespace render
} // namespace ybi
