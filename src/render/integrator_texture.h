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
    const int texelX = ybi::texture::TexelFromUnitUV(uu, textureWidth);
    const int texelY = ybi::texture::TexelFromUnitUV(vv, textureHeight);
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

    if (params.virtualTextureTileEntries == 0ull || params.virtualTextureTilePixels == 0ull ||
        params.virtualTextureTileEntryCapacity <= 0)
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

    const int tileSize = MaxInt(params.feedbackTileSize, 1);
    const int texelX = ybi::texture::TexelFromUnitUV(wrappedU, textureRef.width);
    const int texelY = ybi::texture::TexelFromUnitUV(wrappedV, textureRef.height);
    const unsigned int tileX = ybi::texture::TileCoordFromTexel(texelX, tileSize);
    const unsigned int tileY = ybi::texture::TileCoordFromTexel(texelY, tileSize);
    const int udim = ybi::texture::UdimFromUV(u, v);
    const unsigned int udimBits = ybi::texture::UdimBitsFromUdim(udim);
    const unsigned int textureId =
        static_cast<unsigned int>(ClampInt(geomRef.materialIndex, 0, int((1u << 23u) - 1u)));
    const unsigned int mip = 1u;
    const unsigned long long key =
        ybi::texture::PackVirtualTextureKey(tileX, tileY, udimBits, textureId, mip);

    const LaunchParams::VirtualTextureTileEntry *entries =
        reinterpret_cast<const LaunchParams::VirtualTextureTileEntry *>(
            params.virtualTextureTileEntries);
    const int capacity = params.virtualTextureTileEntryCapacity;
    const unsigned int mask = static_cast<unsigned int>(capacity - 1);
    unsigned int slot = ybi::texture::HashVirtualTextureKey(key) & mask;
    const LaunchParams::VirtualTextureTileEntry *match = nullptr;
    for (int probe = 0; probe < capacity; ++probe)
    {
        const LaunchParams::VirtualTextureTileEntry &entry = entries[slot];
        if (entry.key == key)
        {
            match = &entry;
            break;
        }
        if (entry.key == kVirtualTextureEmptyKey)
        {
            break;
        }
        slot = (slot + 1u) & mask;
    }

    if (!match || match->width == 0u || match->height == 0u)
    {
        outColor = MakeVec3(1.0f, 0.0f, 1.0f);
        return true;
    }

    const int localX = ClampInt(texelX - int(tileX) * tileSize, 0, int(match->width) - 1);
    const int localY = ClampInt(texelY - int(tileY) * tileSize, 0, int(match->height) - 1);
    const unsigned long long sampleOffset =
        match->pixelOffset +
        (static_cast<unsigned long long>(localY) * static_cast<unsigned long long>(match->width) +
         static_cast<unsigned long long>(localX)) *
            4ull;
    const unsigned char *pixels =
        reinterpret_cast<const unsigned char *>(params.virtualTextureTilePixels);
    const Vec4 sample = MakeVec4(float(pixels[sampleOffset + 0]) * (1.0f / 255.0f),
                                 float(pixels[sampleOffset + 1]) * (1.0f / 255.0f),
                                 float(pixels[sampleOffset + 2]) * (1.0f / 255.0f),
                                 float(pixels[sampleOffset + 3]) * (1.0f / 255.0f));
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
