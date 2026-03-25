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

template <typename State>
YBI_DEVICE void TryWriteTextureFeedback(State &state,
                                               const LaunchParams::InstanceGeomRef &geomRef,
                                               int semantic,
                                               unsigned int primitiveIndex,
                                               float wrappedU,
                                               float wrappedV,
                                               int textureWidth,
                                               int textureHeight,
                                               unsigned int udimBits,
                                               unsigned int mip);

template <typename State>
YBI_DEVICE bool TryComputeTextureMipLevelProfiled(State &state,
                                                         const HitInfo &hit,
                                                         int textureWidth,
                                                         int textureHeight,
                                                         unsigned int *outMip)
{
    const unsigned long long start = state.BeginTextureMipTiming();
    const bool result = TryComputeTextureMipLevel(hit, textureWidth, textureHeight, outMip);
    state.EndTextureMipTiming(start, result);
    return result;
}

template <typename State>
YBI_DEVICE void TryWriteTextureFeedbackProfiled(State &state,
                                                       const LaunchParams::InstanceGeomRef &geomRef,
                                                       int semantic,
                                                       unsigned int primitiveIndex,
                                                       float wrappedU,
                                                       float wrappedV,
                                                       int textureWidth,
                                                       int textureHeight,
                                                       unsigned int udimBits,
                                                       unsigned int mip)
{
    const unsigned long long start = state.BeginFeedbackWriteTiming();
    TryWriteTextureFeedback(state,
                            geomRef,
                            semantic,
                            primitiveIndex,
                            wrappedU,
                            wrappedV,
                            textureWidth,
                            textureHeight,
                            udimBits,
                            mip);
    state.EndFeedbackWriteTiming(start);
}

YBI_DEVICE Vec3 MaterialSampleToColorForSemantic(int semantic, const Vec4 &sample)
{
    if (semantic == kSemanticRoughness || semantic == kSemanticMetallic ||
        semantic == kSemanticOcclusion || semantic == kSemanticIor || semantic == kSemanticOpacity ||
        semantic == kSemanticClearcoat || semantic == kSemanticClearcoatRoughness)
    {
        return Vec3(sample.x, sample.x, sample.x);
    }
    return Vec3(sample.x, sample.y, sample.z);
}

YBI_DEVICE Vec3 MaterialSampleToViewColor(const LaunchParams &params, const Vec4 &sample)
{
    return MaterialSampleToColorForSemantic(params.textureViewSemantic, sample);
}

YBI_DEVICE uint32_t PackVirtualTexturePageEntry(unsigned int page,
                                                       unsigned int physicalTextureID,
                                                       unsigned int pageType)
{
    return ((page & 0xffu) << 0u) | ((physicalTextureID & 0x7fffffu) << 8u) | ((pageType & 0x1u) << 31u);
}

YBI_DEVICE void UnpackVirtualTexturePageEntry(uint32_t packed,
                                                     unsigned int *outPage,
                                                     unsigned int *outPhysicalTextureID,
                                                     unsigned int *outPageType)
{
    if (outPage)
    {
        *outPage = (packed >> 0u) & 0xffu;
    }
    if (outPhysicalTextureID)
    {
        *outPhysicalTextureID = (packed >> 8u) & 0x7fffffu;
    }
    if (outPageType)
    {
        *outPageType = (packed >> 31u) & 0x1u;
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

YBI_DEVICE bool IsUsableMaterialTextureRef(const LaunchParams::MaterialTextureRef &ref)
{
    return ref.width > 0 && ref.height > 0;
}

YBI_DEVICE bool ResolveMaterialTextureRefBase(
    const LaunchParams &params,
    const LaunchParams::InstanceGeomRef &geomRef,
    int semantic,
    const LaunchParams::MaterialTextureRef **outMaterialRefs,
    int *outBase,
    int *outMaxSlots)
{
    if (!outMaterialRefs || !outBase || !outMaxSlots || params.materialTextureRefs == 0ull ||
        params.materialTextureRefCount <= 0 || params.materialTextureRefStride <= 0 ||
        geomRef.materialIndex < 0 || geomRef.materialIndex >= params.materialTextureRefCount || semantic < 0 ||
        semantic >= params.materialTextureRefSemanticCount)
    {
        return false;
    }

    *outMaterialRefs =
        reinterpret_cast<const LaunchParams::MaterialTextureRef *>(params.materialTextureRefs);
    *outBase = (geomRef.materialIndex * params.materialTextureRefSemanticCount + semantic) *
               params.materialTextureRefStride;
    *outMaxSlots = params.materialTextureRefCount * params.materialTextureRefSemanticCount *
                   params.materialTextureRefStride;
    return true;
}

YBI_DEVICE bool ResolveMaterialTextureRefBase(
    const LaunchParams &params,
    const LaunchParams::InstanceGeomRef &geomRef,
    const LaunchParams::MaterialTextureRef **outMaterialRefs,
    int *outBase,
    int *outMaxSlots)
{
    return ResolveMaterialTextureRefBase(
        params, geomRef, params.textureViewSemantic, outMaterialRefs, outBase, outMaxSlots);
}

YBI_DEVICE bool FetchMaterialTextureRefForUdimSlot(
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

YBI_DEVICE bool FindFallbackMaterialTextureRef(
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

YBI_DEVICE bool TryComputeMaterialTextureSampleInputs(
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

YBI_DEVICE bool ResolveVirtualTextureTextureMeta(
    const LaunchParams &params,
    const LaunchParams::InstanceGeomRef &geomRef,
    int semantic,
    const LaunchParams::VirtualTextureTextureMeta **outMeta)
{
    if (params.virtualTextureTextureMeta == 0ull || geomRef.materialIndex < 0 || semantic < 0 ||
        semantic >= params.materialTextureRefSemanticCount)
    {
        return false;
    }
    const int textureId = geomRef.materialIndex * params.materialTextureRefSemanticCount + semantic;
    if (textureId < 0 || textureId >= params.virtualTextureTextureMetaCount)
    {
        return false;
    }
    const LaunchParams::VirtualTextureTextureMeta *textures =
        reinterpret_cast<const LaunchParams::VirtualTextureTextureMeta *>(params.virtualTextureTextureMeta);
    if (outMeta)
    {
        *outMeta = &textures[textureId];
    }
    return true;
}

YBI_DEVICE bool TryResolveVirtualTextureLocalUdim(
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

YBI_DEVICE bool ResolveVirtualTextureUdimBits(
    const LaunchParams &params,
    const LaunchParams::InstanceGeomRef &geomRef,
    int semantic,
    float u,
    float v,
    float wrappedU,
    float wrappedV,
    const LaunchParams::VirtualTextureTextureMeta **outMeta,
    unsigned int *outUdimBits,
    unsigned int *outLocalUdim)
{
    const LaunchParams::VirtualTextureTextureMeta *meta = nullptr;
    if (!ResolveVirtualTextureTextureMeta(params, geomRef, semantic, &meta) || !meta)
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

YBI_DEVICE bool TryResolveVirtualTextureTailSample(
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

YBI_DEVICE bool ResolveVirtualTextureInfo(const LaunchParams &params,
                                                 const LaunchParams::InstanceGeomRef &geomRef,
                                                 int semantic,
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
        params.virtualTextureUdimInfos == 0ull)
    {
        return false;
    }

    const LaunchParams::VirtualTextureTextureMeta *metaPtr = nullptr;
    if (!ResolveVirtualTextureTextureMeta(params, geomRef, semantic, &metaPtr) || !metaPtr)
    {
        return false;
    }
    const LaunchParams::VirtualTextureTextureMeta &meta = *metaPtr;
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

YBI_DEVICE bool ReadVirtualTexturePageEntry(const LaunchParams &params,
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
YBI_DEVICE void TryWriteTextureFeedback(State &state,
                                               const LaunchParams::InstanceGeomRef &geomRef,
                                               int semantic,
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
                   374761393u ^
               (static_cast<unsigned int>(ClampInt(semantic, 0, params.materialTextureRefSemanticCount)) +
                1u) *
                   668265263u);
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
    const int textureIdInt = geomRef.materialIndex * params.materialTextureRefSemanticCount + semantic;
    const unsigned int textureId =
        static_cast<unsigned int>(ClampInt(textureIdInt, 0, int((1u << 23u) - 1u)));
    const unsigned long long key =
        ybi::texture::PackVirtualTextureKey(tileX, tileY, udimBits, textureId, mip);

    state.RecordFeedbackKey(key);
}

template <typename State>
YBI_DEVICE bool TrySampleVirtualTexture(State &state,
                                               const LaunchParams::InstanceGeomRef &geomRef,
                                               int semantic,
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

    if (params.virtualTexturePageTableEntries == 0ull || params.virtualTextureTextureMeta == 0ull ||
        params.virtualTextureMipInfos == 0ull)
    {
        outColor = Vec3(0.0f, 0.0f, 0.0f);
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
    if (!ResolveVirtualTextureTextureMeta(params, geomRef, semantic, &meta) ||
        !TryResolveVirtualTextureLocalUdim(*meta, udimBits, &localUdim))
    {
        outColor = Vec3(0.0f, 0.0f, 0.0f);
        return true;
    }

    const LaunchParams::VirtualTextureMipInfo *mipInfo = nullptr;
    unsigned int vaX = 0u;
    unsigned int vaY = 0u;
    Vec4 sample = Vec4(0.0f, 0.0f, 0.0f, 0.0f);
    bool haveSample = false;
    const unsigned char *samplePixels = nullptr;
    unsigned long long sampleOffset = 0ull;
    if (!ResolveVirtualTextureInfo(
            params, geomRef, semantic, mip, udimBits, tileX, tileY, &meta, &mipInfo, &localUdim, &vaX, &vaY))
    {
        if (!TryResolveVirtualTextureTailSample(
                meta, localUdim, wrappedU, wrappedV, tileSize, &samplePixels, &sampleOffset))
        {
            outColor = Vec3(0.0f, 0.0f, 0.0f);
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
                outColor = Vec3(0.0f, 0.0f, 0.0f);
                return true;
            }
        }
        else
        {
            unsigned int page = 0u;
            unsigned int physicalTextureID = 0u;
            unsigned int pageType = 0u;
            UnpackVirtualTexturePageEntry(packedEntry, &page, &physicalTextureID, &pageType);

            const int localX = ClampInt(texelX - int(tileX) * tileSize, 0, tileSize - 1);
            const int localY = ClampInt(texelY - int(tileY) * tileSize, 0, tileSize - 1);
            if (pageType == kVirtualTexturePageTypeStream &&
                page < static_cast<unsigned int>(MaxInt(params.virtualTexturePhysicalPagesPerTexture, 0)) &&
                physicalTextureID <
                    static_cast<unsigned int>(MaxInt(params.virtualTexturePhysicalTextureCount, 0)))
            {
                haveSample = state.SampleVirtualTexturePage(physicalTextureID, page, localX, localY, sample);
            }
            else if (pageType == kVirtualTexturePageTypeTail && meta->tailPixels != 0ull &&
                     page < meta->tailPageCountX)
            {
                const int tailX = ybi::texture::TexelFromUnitUV(wrappedU, tileSize);
                const int tailY = ybi::texture::TexelFromUnitUV(wrappedV, tileSize);
                const unsigned long long pageIndex =
                    static_cast<unsigned long long>(page);
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
        if (!haveSample)
        {
            outColor = Vec3(0.0f, 0.0f, 0.0f);
            return true;
        }
    }

    if (!haveSample)
    {
        sample = Vec4(float(samplePixels[sampleOffset + 0]) * (1.0f / 255.0f),
                      float(samplePixels[sampleOffset + 1]) * (1.0f / 255.0f),
                      float(samplePixels[sampleOffset + 2]) * (1.0f / 255.0f),
                      float(samplePixels[sampleOffset + 3]) * (1.0f / 255.0f));
    }
    outColor = MaterialSampleToViewColor(params, sample);
    return true;
}

template <typename State>
YBI_DEVICE bool TrySampleMaterialTexture(State &state,
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
        outColor = Vec3(0.0f, 0.0f, 0.0f);
        return true;
    }

    if (params.virtualTextureEnabled != 0)
    {
        unsigned int resolvedUdimBits = 0u;
        if (!ResolveVirtualTextureUdimBits(params,
                                           geomRef,
                                           params.textureViewSemantic,
                                           inputs.u,
                                           inputs.v,
                                           wrappedU,
                                           wrappedV,
                                           nullptr,
                                           &resolvedUdimBits,
                                           nullptr))
        {
            outColor = Vec3(0.0f, 0.0f, 0.0f);
            return true;
        }

        LaunchParams::MaterialTextureRef vtRef = {};
        const int resolvedUdimSlot = ybi::texture::UdimSlotFromUdim(
            1001 + int(resolvedUdimBits), params.materialTextureRefStride);
        if (!FetchMaterialTextureRefForUdimSlot(materialRefs, base, maxSlots, resolvedUdimSlot, &vtRef) ||
            !IsUsableMaterialTextureRef(vtRef))
        {
            outColor = Vec3(0.0f, 0.0f, 0.0f);
            return true;
        }

        unsigned int sampleMip = 0u;
        const bool haveMip =
            TryComputeTextureMipLevelProfiled(state, hit, vtRef.width, vtRef.height, &sampleMip);
        if (!haveMip)
        {
            sampleMip = static_cast<unsigned int>(
                ClampInt(MaxInt(params.virtualTextureSampleMip, 0),
                         0,
                         ComputeTextureMipCount(vtRef.width, vtRef.height) - 1));
        }

        if (TrySampleVirtualTexture(
                state,
                geomRef,
                params.textureViewSemantic,
                vtRef,
                wrappedU,
                wrappedV,
                resolvedUdimBits,
                sampleMip,
                outColor))
        {
            if (haveMip)
            {
                TryWriteTextureFeedbackProfiled(state,
                                                geomRef,
                                                params.textureViewSemantic,
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
    if (TryComputeTextureMipLevelProfiled(state, hit, rawRef.width, rawRef.height, &feedbackMip))
    {
        TryWriteTextureFeedbackProfiled(state,
                                        geomRef,
                                        params.textureViewSemantic,
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

template <typename State>
YBI_DEVICE bool TrySampleMaterialTextureSemantic(State &state,
                                                        const LaunchParams::InstanceGeomRef &geomRef,
                                                        const HitInfo &hit,
                                                        int semantic,
                                                        Vec4 &outSample)
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
    if (!ResolveMaterialTextureRefBase(params, geomRef, semantic, &materialRefs, &base, &maxSlots))
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
        outSample = Vec4(0.0f, 0.0f, 0.0f, 0.0f);
        return true;
    }

    if (params.virtualTextureEnabled != 0)
    {
        Vec3 vtColor = {};
        unsigned int resolvedUdimBits = 0u;
        if (ResolveVirtualTextureUdimBits(params,
                                          geomRef,
                                          semantic,
                                          inputs.u,
                                          inputs.v,
                                          wrappedU,
                                          wrappedV,
                                          nullptr,
                                          &resolvedUdimBits,
                                          nullptr))
        {
            LaunchParams::MaterialTextureRef vtRef = {};
            const int resolvedUdimSlot = ybi::texture::UdimSlotFromUdim(
                1001 + int(resolvedUdimBits), params.materialTextureRefStride);
            if (FetchMaterialTextureRefForUdimSlot(materialRefs, base, maxSlots, resolvedUdimSlot, &vtRef) &&
                IsUsableMaterialTextureRef(vtRef))
            {
                unsigned int sampleMip = 0u;
                const bool haveMip = TryComputeTextureMipLevelProfiled(
                    state, hit, vtRef.width, vtRef.height, &sampleMip);
                if (!haveMip)
                {
                    sampleMip = static_cast<unsigned int>(
                        ClampInt(MaxInt(params.virtualTextureSampleMip, 0),
                                 0,
                                 ComputeTextureMipCount(vtRef.width, vtRef.height) - 1));
                }

                if (TrySampleVirtualTexture(
                        state,
                        geomRef,
                        semantic,
                        vtRef,
                        wrappedU,
                        wrappedV,
                        resolvedUdimBits,
                        sampleMip,
                        vtColor))
                {
                    outSample = Vec4(vtColor.x, vtColor.y, vtColor.z, 1.0f);
                    if (haveMip)
                    {
                        TryWriteTextureFeedbackProfiled(state,
                                                        geomRef,
                                                        semantic,
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

    if (!state.SampleTexture2D(rawRef, inputs.unitU, inputs.unitV, outSample))
    {
        return false;
    }
    unsigned int feedbackMip = 0u;
    if (TryComputeTextureMipLevelProfiled(state, hit, rawRef.width, rawRef.height, &feedbackMip))
    {
        TryWriteTextureFeedbackProfiled(state,
                                        geomRef,
                                        semantic,
                                        static_cast<unsigned int>(hit.primitiveIndex),
                                        inputs.unitU,
                                        inputs.unitV,
                                        rawRef.width,
                                        rawRef.height,
                                        inputs.rawUdimBits,
                                        feedbackMip);
    }
    return true;
}

} // namespace integrator
} // namespace render
} // namespace ybi
