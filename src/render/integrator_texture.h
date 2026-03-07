#pragma once

#include "render/integrator_common.h"
#include "render/integrator_texture_mip.h"
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

struct MaterialTextureSampleInputs
{
    UV2 uv0 = {0.0f, 0.0f};
    UV2 uv1 = {0.0f, 0.0f};
    UV2 uv2 = {0.0f, 0.0f};
    float u = 0.0f;
    float v = 0.0f;
    float unitU = 0.0f;
    float unitV = 0.0f;
    int rawUdim = 1001;
    unsigned int rawUdimBits = 0u;
};

YBI_INTEGRATOR_HD bool IsUsableMaterialTextureRef(const LaunchParams::MaterialTextureRef &ref)
{
    return ref.width > 0 && ref.height > 0;
}

YBI_INTEGRATOR_HD bool ResolveMaterialTextureRefBase(
    const LaunchParams &params,
    const LaunchParams::InstanceGeomRef &geomRef,
    const LaunchParams::MaterialTextureRef **outMaterialRefs,
    int *outBase,
    int *outMaxSlots)
{
    if (!outMaterialRefs || !outBase || !outMaxSlots || params.materialTextureRefs == 0ull ||
        params.materialTextureRefCount <= 0 || params.materialTextureRefStride <= 0 ||
        geomRef.materialIndex < 0 || geomRef.materialIndex >= params.materialTextureRefCount ||
        params.textureViewSemantic < 0 ||
        params.textureViewSemantic >= params.materialTextureRefSemanticCount)
    {
        return false;
    }

    *outMaterialRefs =
        reinterpret_cast<const LaunchParams::MaterialTextureRef *>(params.materialTextureRefs);
    *outBase = (geomRef.materialIndex * params.materialTextureRefSemanticCount +
                params.textureViewSemantic) *
               params.materialTextureRefStride;
    *outMaxSlots = params.materialTextureRefCount * params.materialTextureRefSemanticCount *
                   params.materialTextureRefStride;
    return true;
}

YBI_INTEGRATOR_HD bool FetchMaterialTextureRefForUdimSlot(
    const LaunchParams::MaterialTextureRef *materialRefs,
    int base,
    int maxSlots,
    int udimSlot,
    LaunchParams::MaterialTextureRef *outRef)
{
    if (!materialRefs || !outRef)
    {
        return false;
    }

    const int slot = base + udimSlot;
    if (slot < 0 || slot >= maxSlots)
    {
        return false;
    }

    *outRef = materialRefs[slot];
    return true;
}

YBI_INTEGRATOR_HD bool FindFallbackMaterialTextureRef(
    const LaunchParams &params,
    const LaunchParams::MaterialTextureRef *materialRefs,
    int base,
    int maxSlots,
    LaunchParams::MaterialTextureRef *outRef)
{
    if (!materialRefs || !outRef)
    {
        return false;
    }

    for (int udimSlot = 0; udimSlot < params.materialTextureRefStride; ++udimSlot)
    {
        LaunchParams::MaterialTextureRef candidate = {};
        if (!FetchMaterialTextureRefForUdimSlot(materialRefs, base, maxSlots, udimSlot, &candidate))
        {
            continue;
        }
        if (IsUsableMaterialTextureRef(candidate))
        {
            *outRef = candidate;
            return true;
        }
    }
    return false;
}

YBI_INTEGRATOR_HD bool TryComputeMaterialTextureSampleInputs(
    const LaunchParams::InstanceGeomRef &geomRef,
    unsigned int primitiveIndex,
    const Vec3 &barycentrics,
    MaterialTextureSampleInputs *outInputs)
{
    if (!outInputs || geomRef.texcoords == 0ull || geomRef.texcoordIndices == 0ull)
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
    outInputs->uv0 = texcoords[t0];
    outInputs->uv1 = texcoords[t1];
    outInputs->uv2 = texcoords[t2];
    outInputs->u = outInputs->uv0.x * barycentrics.x + outInputs->uv1.x * barycentrics.y +
                   outInputs->uv2.x * barycentrics.z;
    outInputs->v = outInputs->uv0.y * barycentrics.x + outInputs->uv1.y * barycentrics.y +
                   outInputs->uv2.y * barycentrics.z;
    outInputs->unitU = ybi::texture::UnitUV(outInputs->u);
    outInputs->unitV = ybi::texture::UnitUV(outInputs->v);
    outInputs->rawUdim = ybi::texture::UdimFromUV(outInputs->u, outInputs->v);
    outInputs->rawUdimBits = ybi::texture::UdimBitsFromUdim(outInputs->rawUdim);
    return true;
}

YBI_INTEGRATOR_HD bool ResolveVirtualTextureTextureMeta(
    const LaunchParams &params,
    const LaunchParams::InstanceGeomRef &geomRef,
    const LaunchParams::VirtualTextureTextureMeta **outMeta)
{
    if (params.virtualTextureTextureMeta == 0ull || geomRef.materialIndex < 0 ||
        geomRef.materialIndex >= params.virtualTextureTextureMetaCount)
    {
        return false;
    }
    const LaunchParams::VirtualTextureTextureMeta *textures =
        reinterpret_cast<const LaunchParams::VirtualTextureTextureMeta *>(params.virtualTextureTextureMeta);
    if (outMeta)
    {
        *outMeta = &textures[geomRef.materialIndex];
    }
    return true;
}

YBI_INTEGRATOR_HD bool TryResolveVirtualTextureLocalUdim(
    const LaunchParams::VirtualTextureTextureMeta &meta,
    unsigned int udimBits,
    unsigned int *outLocalUdim)
{
    if (udimBits >= 128u)
    {
        return false;
    }
    const int local = static_cast<int>(meta.udimToLocal[udimBits]);
    if (local < 0)
    {
        return false;
    }
    if (outLocalUdim)
    {
        *outLocalUdim = static_cast<unsigned int>(local);
    }
    return true;
}

YBI_INTEGRATOR_HD bool ResolveVirtualTextureUdimBits(
    const LaunchParams &params,
    const LaunchParams::InstanceGeomRef &geomRef,
    float u,
    float v,
    float wrappedU,
    float wrappedV,
    const LaunchParams::VirtualTextureTextureMeta **outMeta,
    unsigned int *outUdimBits,
    unsigned int *outLocalUdim)
{
    const LaunchParams::VirtualTextureTextureMeta *meta = nullptr;
    if (!ResolveVirtualTextureTextureMeta(params, geomRef, &meta) || !meta)
    {
        return false;
    }

    const unsigned int rawBits =
        ybi::texture::UdimBitsFromUdim(ybi::texture::UdimFromUV(u, v));
    unsigned int localUdim = 0u;
    if (TryResolveVirtualTextureLocalUdim(*meta, rawBits, &localUdim))
    {
        if (outMeta)
        {
            *outMeta = meta;
        }
        if (outUdimBits)
        {
            *outUdimBits = rawBits;
        }
        if (outLocalUdim)
        {
            *outLocalUdim = localUdim;
        }
        return true;
    }

    const unsigned int wrappedBits =
        ybi::texture::UdimBitsFromUdim(ybi::texture::UdimFromUV(wrappedU, wrappedV));
    if (TryResolveVirtualTextureLocalUdim(*meta, wrappedBits, &localUdim))
    {
        if (outMeta)
        {
            *outMeta = meta;
        }
        if (outUdimBits)
        {
            *outUdimBits = wrappedBits;
        }
        if (outLocalUdim)
        {
            *outLocalUdim = localUdim;
        }
        return true;
    }
    return false;
}

YBI_INTEGRATOR_HD bool TryResolveVirtualTextureTailSample(
    const LaunchParams::VirtualTextureTextureMeta *meta,
    unsigned int localUdim,
    float wrappedU,
    float wrappedV,
    int tileSize,
    const unsigned char **outSamplePixels,
    unsigned long long *outSampleOffset)
{
    if (!meta || !outSamplePixels || !outSampleOffset || meta->tailPixels == 0ull || meta->tailPageCountX == 0u ||
        meta->tailPageCountY == 0u)
    {
        return false;
    }

    const unsigned long long pageCapacity = static_cast<unsigned long long>(meta->tailPageCountX) *
                                            static_cast<unsigned long long>(meta->tailPageCountY);
    if (static_cast<unsigned long long>(localUdim) >= pageCapacity)
    {
        return false;
    }

    const int safeTileSize = MaxInt(tileSize, 1);
    const int tailX = ybi::texture::TexelFromUnitUV(wrappedU, safeTileSize);
    const int tailY = ybi::texture::TexelFromUnitUV(wrappedV, safeTileSize);
    const unsigned int fallbackPageX = localUdim % meta->tailPageCountX;
    const unsigned int fallbackPageY = localUdim / meta->tailPageCountX;
    if (fallbackPageY >= meta->tailPageCountY)
    {
        return false;
    }

    const unsigned long long pageIndex =
        static_cast<unsigned long long>(fallbackPageY) * static_cast<unsigned long long>(meta->tailPageCountX) +
        static_cast<unsigned long long>(fallbackPageX);
    const unsigned long long pageBytes =
        static_cast<unsigned long long>(safeTileSize) * static_cast<unsigned long long>(safeTileSize) * 4ull;
    const unsigned long long sampleOffset =
        pageIndex * pageBytes +
        (static_cast<unsigned long long>(tailY) * static_cast<unsigned long long>(safeTileSize) +
         static_cast<unsigned long long>(tailX)) *
            4ull;
    *outSamplePixels = reinterpret_cast<const unsigned char *>(meta->tailPixels);
    *outSampleOffset = sampleOffset;
    return true;
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
        params.virtualTextureUdimInfos == 0ull ||
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
    const unsigned int mipInfoIndex = meta.mipInfoOffset + mip;
    if (mipInfoIndex >= static_cast<unsigned int>(MaxInt(params.virtualTextureMipInfoCount, 0)))
    {
        return false;
    }
    const LaunchParams::VirtualTextureMipInfo &mipInfo = mipInfos[mipInfoIndex];
    if (static_cast<unsigned int>(localUdim) >= mipInfo.udimInfoCount)
    {
        return false;
    }
    const LaunchParams::VirtualTextureUdimInfo *udimInfos =
        reinterpret_cast<const LaunchParams::VirtualTextureUdimInfo *>(params.virtualTextureUdimInfos);
    const unsigned int udimInfoIndex = mipInfo.udimInfoOffset + static_cast<unsigned int>(localUdim);
    if (udimInfoIndex >= static_cast<unsigned int>(MaxInt(params.virtualTextureUdimInfoCount, 0)))
    {
        return false;
    }
    const LaunchParams::VirtualTextureUdimInfo &udimInfo =
        udimInfos[udimInfoIndex];
    if (tileX >= udimInfo.pageCountX || tileY >= udimInfo.pageCountY)
    {
        return false;
    }

    const unsigned int vaX = udimInfo.basePageX + tileX;
    const unsigned int vaY = udimInfo.basePageY + tileY;
    if (vaX >= udimInfo.basePageX + udimInfo.pageCountX || vaY >= udimInfo.basePageY + udimInfo.pageCountY)
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
                                               float wrappedU,
                                               float wrappedV,
                                               int textureWidth,
                                               int textureHeight,
                                               unsigned int udimBits,
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
    const int texelX = ybi::texture::TexelFromUnitUV(wrappedU, mipWidth);
    const int texelY = ybi::texture::TexelFromUnitUV(wrappedV, mipHeight);
    const unsigned int tileX = ybi::texture::TileCoordFromTexel(texelX, tileSize);
    const unsigned int tileY = ybi::texture::TileCoordFromTexel(texelY, tileSize);
    const unsigned int textureId =
        static_cast<unsigned int>(ClampInt(geomRef.materialIndex, 0, int((1u << 23u) - 1u)));
    const unsigned long long key =
        ybi::texture::PackVirtualTextureKey(tileX, tileY, udimBits, textureId, mip);

    state.RecordFeedbackKey(key);
}

template <typename State>
YBI_INTEGRATOR_HD bool TrySampleVirtualTexture(State &state,
                                               const LaunchParams::InstanceGeomRef &geomRef,
                                               const LaunchParams::MaterialTextureRef &textureRef,
                                               float wrappedU,
                                               float wrappedV,
                                               unsigned int udimBits,
                                               unsigned int mip,
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
        outColor = MakeVec3(0.0f, 0.0f, 0.0f);
        return true;
    }

    const int mipWidth = MaxInt(textureRef.width >> int(mip), 1);
    const int mipHeight = MaxInt(textureRef.height >> int(mip), 1);
    const int tileSize = MaxInt(params.virtualTexturePageSize, 1);
    const int texelX = ybi::texture::TexelFromUnitUV(wrappedU, mipWidth);
    const int texelY = ybi::texture::TexelFromUnitUV(wrappedV, mipHeight);
    const unsigned int tileX = ybi::texture::TileCoordFromTexel(texelX, tileSize);
    const unsigned int tileY = ybi::texture::TileCoordFromTexel(texelY, tileSize);
    const LaunchParams::VirtualTextureTextureMeta *meta = nullptr;
    unsigned int localUdim = 0u;
    if (!ResolveVirtualTextureTextureMeta(params, geomRef, &meta) ||
        !TryResolveVirtualTextureLocalUdim(*meta, udimBits, &localUdim))
    {
        outColor = MakeVec3(0.0f, 0.0f, 0.0f);
        return true;
    }

    const LaunchParams::VirtualTextureMipInfo *mipInfo = nullptr;
    unsigned int vaX = 0u;
    unsigned int vaY = 0u;
    const unsigned char *samplePixels = nullptr;
    unsigned long long sampleOffset = 0ull;
    if (!ResolveVirtualTextureInfo(
            params, geomRef, mip, udimBits, tileX, tileY, &meta, &mipInfo, &localUdim, &vaX, &vaY))
    {
        if (!TryResolveVirtualTextureTailSample(
                meta, localUdim, wrappedU, wrappedV, tileSize, &samplePixels, &sampleOffset))
        {
            outColor = MakeVec3(0.0f, 0.0f, 0.0f);
            return true;
        }
    }
    else
    {
        uint32_t packedEntry = 0u;
        if (!ReadVirtualTexturePageEntry(params, mipInfo->level, vaX, vaY, &packedEntry))
        {
            if (!TryResolveVirtualTextureTailSample(
                    meta, localUdim, wrappedU, wrappedV, tileSize, &samplePixels, &sampleOffset))
            {
                outColor = MakeVec3(0.0f, 0.0f, 0.0f);
                return true;
            }
        }
        else
        {
            unsigned int pageX = 0u;
            unsigned int pageY = 0u;
            unsigned int pageType = 0u;
            unsigned int flags = 0u;
            UnpackVirtualTexturePageEntry(packedEntry, &pageX, &pageY, &pageType, &flags);

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
            else
            {
                TryResolveVirtualTextureTailSample(
                    meta, localUdim, wrappedU, wrappedV, tileSize, &samplePixels, &sampleOffset);
            }
        }
    }

    if (!samplePixels)
    {
        outColor = MakeVec3(0.0f, 0.0f, 0.0f);
        return true;
    }

    const Vec4 sample = MakeVec4(float(samplePixels[sampleOffset + 0]) * (1.0f / 255.0f),
                                 float(samplePixels[sampleOffset + 1]) * (1.0f / 255.0f),
                                 float(samplePixels[sampleOffset + 2]) * (1.0f / 255.0f),
                                 float(samplePixels[sampleOffset + 3]) * (1.0f / 255.0f));
    outColor = MaterialSampleToViewColor(params, sample);
    return true;
}

template <typename State>
YBI_INTEGRATOR_HD bool
TryWriteFeedbackOnly(State &state,
                     const LaunchParams::InstanceGeomRef &geomRef,
                     const HitInfo &hit)
{
    const LaunchParams &params = state.Params();
    if (params.materialTextureRefs == 0ull || params.materialTextureRefCount <= 0 ||
        params.materialTextureRefStride <= 0 || !hit.hasBarycentrics)
    {
        return false;
    }

    MaterialTextureSampleInputs inputs = {};
    if (!TryComputeMaterialTextureSampleInputs(
            geomRef, static_cast<unsigned int>(hit.primitiveIndex), hit.barycentrics, &inputs))
    {
        return false;
    }

    const LaunchParams::MaterialTextureRef *materialRefs = nullptr;
    int base = 0;
    int maxSlots = 0;
    if (!ResolveMaterialTextureRefBase(params, geomRef, &materialRefs, &base, &maxSlots))
    {
        return false;
    }

    const int rawUdimSlot =
        ybi::texture::UdimSlotFromUdim(inputs.rawUdim, params.materialTextureRefStride);
    LaunchParams::MaterialTextureRef rawRef = {};
    const bool hasRawRef =
        FetchMaterialTextureRefForUdimSlot(materialRefs, base, maxSlots, rawUdimSlot, &rawRef);

    LaunchParams::MaterialTextureRef wrapRef = rawRef;
    if (!hasRawRef || !IsUsableMaterialTextureRef(wrapRef))
    {
        if (!FindFallbackMaterialTextureRef(params, materialRefs, base, maxSlots, &wrapRef))
        {
            return false;
        }
    }

    bool blackS = false;
    bool blackT = false;
    const float wrappedU = ApplyWrapMode(inputs.unitU, wrapRef.wrapS, blackS);
    const float wrappedV = ApplyWrapMode(inputs.unitV, wrapRef.wrapT, blackT);
    if (blackS || blackT)
    {
        return false;
    }

    unsigned int feedbackUdimBits = inputs.rawUdimBits;
    LaunchParams::MaterialTextureRef feedbackRef = rawRef;
    if (params.virtualTextureEnabled != 0)
    {
        unsigned int resolvedUdimBits = 0u;
        if (!ResolveVirtualTextureUdimBits(params,
                                           geomRef,
                                           inputs.u,
                                           inputs.v,
                                           wrappedU,
                                           wrappedV,
                                           nullptr,
                                           &resolvedUdimBits,
                                           nullptr))
        {
            return false;
        }

        const int resolvedUdimSlot = ybi::texture::UdimSlotFromUdim(
            1001 + int(resolvedUdimBits), params.materialTextureRefStride);
        if (!FetchMaterialTextureRefForUdimSlot(
                materialRefs, base, maxSlots, resolvedUdimSlot, &feedbackRef) ||
            !IsUsableMaterialTextureRef(feedbackRef))
        {
            return false;
        }
        feedbackUdimBits = resolvedUdimBits;
    }
    else if (!hasRawRef || !IsUsableMaterialTextureRef(feedbackRef))
    {
        return false;
    }

    unsigned int mip = 0u;
    if (!TryComputeTextureMipLevel(
            params, hit, inputs.uv0, inputs.uv1, inputs.uv2, feedbackRef.width, feedbackRef.height, &mip))
    {
        return false;
    }

    TryWriteTextureFeedback(
        state,
        geomRef,
        static_cast<unsigned int>(hit.primitiveIndex),
        wrappedU,
        wrappedV,
        feedbackRef.width,
        feedbackRef.height,
        feedbackUdimBits,
        mip);
    return true;
}

template <typename State>
YBI_INTEGRATOR_HD bool TrySampleMaterialTexture(State &state,
                                                const LaunchParams::InstanceGeomRef &geomRef,
                                                const HitInfo &hit,
                                                Vec3 &outColor)
{
    const LaunchParams &params = state.Params();
    if (params.materialTextureRefs == 0ull || params.materialTextureRefCount <= 0 ||
        params.materialTextureRefStride <= 0 || !hit.hasBarycentrics)
    {
        return false;
    }

    MaterialTextureSampleInputs inputs = {};
    if (!TryComputeMaterialTextureSampleInputs(
            geomRef, static_cast<unsigned int>(hit.primitiveIndex), hit.barycentrics, &inputs))
    {
        return false;
    }

    const LaunchParams::MaterialTextureRef *materialRefs = nullptr;
    int base = 0;
    int maxSlots = 0;
    if (!ResolveMaterialTextureRefBase(params, geomRef, &materialRefs, &base, &maxSlots))
    {
        return false;
    }

    const int rawUdimSlot =
        ybi::texture::UdimSlotFromUdim(inputs.rawUdim, params.materialTextureRefStride);
    LaunchParams::MaterialTextureRef rawRef = {};
    const bool hasRawRef =
        FetchMaterialTextureRefForUdimSlot(materialRefs, base, maxSlots, rawUdimSlot, &rawRef);
    LaunchParams::MaterialTextureRef wrapRef = rawRef;
    if (!hasRawRef || !IsUsableMaterialTextureRef(wrapRef))
    {
        if (!FindFallbackMaterialTextureRef(params, materialRefs, base, maxSlots, &wrapRef))
        {
            return false;
        }
    }

    bool blackS = false;
    bool blackT = false;
    const float wrappedU = ApplyWrapMode(inputs.unitU, wrapRef.wrapS, blackS);
    const float wrappedV = ApplyWrapMode(inputs.unitV, wrapRef.wrapT, blackT);
    if (blackS || blackT)
    {
        outColor = MakeVec3(0.0f, 0.0f, 0.0f);
        return true;
    }

    if (params.virtualTextureEnabled != 0)
    {
        unsigned int resolvedUdimBits = 0u;
        if (!ResolveVirtualTextureUdimBits(params,
                                           geomRef,
                                           inputs.u,
                                           inputs.v,
                                           wrappedU,
                                           wrappedV,
                                           nullptr,
                                           &resolvedUdimBits,
                                           nullptr))
        {
            outColor = MakeVec3(0.0f, 0.0f, 0.0f);
            return true;
        }

        LaunchParams::MaterialTextureRef vtRef = {};
        const int resolvedUdimSlot = ybi::texture::UdimSlotFromUdim(
            1001 + int(resolvedUdimBits), params.materialTextureRefStride);
        if (!FetchMaterialTextureRefForUdimSlot(materialRefs, base, maxSlots, resolvedUdimSlot, &vtRef) ||
            !IsUsableMaterialTextureRef(vtRef))
        {
            outColor = MakeVec3(0.0f, 0.0f, 0.0f);
            return true;
        }

        unsigned int sampleMip = 0u;
        const bool haveMip = TryComputeTextureMipLevel(
            params, hit, inputs.uv0, inputs.uv1, inputs.uv2, vtRef.width, vtRef.height, &sampleMip);
        if (!haveMip)
        {
            sampleMip = static_cast<unsigned int>(
                ClampInt(MaxInt(params.virtualTextureSampleMip, 0),
                         0,
                         ComputeTextureMipCount(vtRef.width, vtRef.height) - 1));
        }

        if (TrySampleVirtualTexture(
                state, geomRef, vtRef, wrappedU, wrappedV, resolvedUdimBits, sampleMip, outColor))
        {
            if (haveMip)
            {
                TryWriteTextureFeedback(state,
                                        geomRef,
                                        static_cast<unsigned int>(hit.primitiveIndex),
                                        wrappedU,
                                        wrappedV,
                                        vtRef.width,
                                        vtRef.height,
                                        resolvedUdimBits,
                                        sampleMip);
            }
            return true;
        }
    }

    if (!hasRawRef || !IsUsableMaterialTextureRef(rawRef))
    {
        return false;
    }

    if (rawRef.textureObject == 0ull || rawRef.valid == 0)
    {
        return false;
    }

    Vec4 sample = {};
    if (!state.SampleTexture2D(rawRef, inputs.unitU, inputs.unitV, sample))
    {
        return false;
    }
    outColor = MaterialSampleToViewColor(params, sample);
    unsigned int feedbackMip = 0u;
    if (TryComputeTextureMipLevel(
            params, hit, inputs.uv0, inputs.uv1, inputs.uv2, rawRef.width, rawRef.height, &feedbackMip))
    {
        TryWriteTextureFeedback(state,
                                geomRef,
                                static_cast<unsigned int>(hit.primitiveIndex),
                                inputs.unitU,
                                inputs.unitV,
                                rawRef.width,
                                rawRef.height,
                                inputs.rawUdimBits,
                                feedbackMip);
    }
    state.MaybeLogSampleSuccess();
    return true;
}

} // namespace integrator
} // namespace render
} // namespace ybi
