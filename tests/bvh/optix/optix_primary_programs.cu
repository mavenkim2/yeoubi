#include <optix.h>
#include <optix_device.h>
#include <assert.h>
#include "texture/virtual_texture_key.h"

struct LaunchParams
{
    struct InstanceGeomRef
    {
        unsigned long long positions;
        unsigned long long indices;
        unsigned long long texcoords;
        unsigned long long texcoordIndices;
        int numPositions;
        int numIndices;
        int numTexcoords;
        int numTexcoordIndices;
        int materialIndex;
    };

    struct WireframeConfig
    {
        float lineWidth;
        float lineFeather;
        float edgeDarkness;
        float padding;
    };

    struct MaterialTextureRef
    {
        unsigned long long textureObject;
        int width;
        int height;
        int valid;
        int wrapS;
        int wrapT;
        int _padding0;
        int _padding1;
    };

    struct VirtualTextureTileEntry
    {
        unsigned long long key;
        unsigned long long pixelOffset;
        unsigned int width;
        unsigned int height;
    };

    OptixTraversableHandle traversable;
    unsigned long long image;
    int width;
    int height;
    float3 cameraOrigin;
    float3 cameraU;
    float3 cameraV;
    float3 cameraW;
    WireframeConfig wireframe;
    int integrator;
    int spp;
    float aoBias;
    float aoMaxDistance;
    unsigned long long instanceGeomRefs;
    int instanceGeomRefCount;
    unsigned long long materialTextureRefs;
    int materialTextureRefCount;
    int materialTextureRefStride;
    int materialTextureRefSemanticCount;
    int textureViewSemantic;
    unsigned long long feedbackKeys;
    unsigned long long feedbackStats;
    int feedbackCapacity;
    int feedbackSamplePercent;
    int feedbackTileSize;
    int currentSpp;
    unsigned long long virtualTextureTileEntries;
    unsigned long long virtualTextureTilePixels;
    int virtualTextureTileEntryCapacity;
    int virtualTextureEnabled;
};

struct HitgroupData
{
    unsigned long long positions;
    unsigned long long indices;
    int numPositions;
    int numIndices;
};

struct float2_simple
{
    float x;
    float y;
};

extern "C"
{
    __constant__ LaunchParams params;
}

static constexpr int kSemanticRoughness = 1;
static constexpr int kSemanticMetallic = 2;
static constexpr int kSemanticOcclusion = 3;
static constexpr int kSemanticIor = 5;
static constexpr int kSemanticOpacity = 7;
static constexpr int kWrapModeUnknown = 0;
static constexpr int kWrapModeRepeat = 1;
static constexpr int kWrapModeClamp = 2;
static constexpr int kWrapModeMirror = 3;
static constexpr int kWrapModeBlack = 4;
static constexpr int kWrapModeUseMetadata = 5;
static constexpr unsigned long long kVirtualTextureEmptyKey = ~0ull;

__device__ unsigned int g_try_sample_logged = 0u;
__device__ unsigned int g_try_sample_fail_mask = 0u;

static __forceinline__ __device__ float3 Normalize3(const float3 &v)
{
    const float invLen = rsqrtf(v.x * v.x + v.y * v.y + v.z * v.z + 1e-20f);
    return make_float3(v.x * invLen, v.y * invLen, v.z * invLen);
}

static __forceinline__ __device__ float3 Cross3(const float3 &a, const float3 &b)
{
    return make_float3(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
}

static __forceinline__ __device__ unsigned int Hash32(unsigned int x)
{
    x ^= x >> 16;
    x *= 0x7feb352du;
    x ^= x >> 15;
    x *= 0x846ca68bu;
    x ^= x >> 16;
    return x;
}

static __forceinline__ __device__ float Random01(unsigned int &state)
{
    state = Hash32(state + 0x9e3779b9u);
    return float(state & 0x00ffffffu) / float(0x01000000u);
}

static __forceinline__ __device__ void BuildOrthonormalBasis(const float3 &n, float3 &t, float3 &b)
{
    const float3 up =
        fabsf(n.z) < 0.999f ? make_float3(0.0f, 0.0f, 1.0f) : make_float3(0.0f, 1.0f, 0.0f);
    t = Normalize3(Cross3(up, n));
    b = Normalize3(Cross3(n, t));
}

static __forceinline__ __device__ float3 SampleCosineHemisphere(float u1, float u2)
{
    const float r = sqrtf(fmaxf(u1, 0.0f));
    const float phi = 6.28318530718f * u2;
    const float x = r * cosf(phi);
    const float y = r * sinf(phi);
    const float z = sqrtf(fmaxf(1.0f - u1, 0.0f));
    return make_float3(x, y, z);
}

static __forceinline__ __device__ float3 FaceForward(const float3 &normal,
                                                     const float3 &referenceDirection)
{
    const float d = normal.x * referenceDirection.x + normal.y * referenceDirection.y +
                    normal.z * referenceDirection.z;
    if (d < 0.0f)
    {
        return normal;
    }
    return make_float3(-normal.x, -normal.y, -normal.z);
}

static __forceinline__ __device__ unsigned int PackColor(float r, float g, float b)
{
    const unsigned int ru = (unsigned int)(fminf(fmaxf(r, 0.0f), 1.0f) * 255.0f + 0.5f);
    const unsigned int gu = (unsigned int)(fminf(fmaxf(g, 0.0f), 1.0f) * 255.0f + 0.5f);
    const unsigned int bu = (unsigned int)(fminf(fmaxf(b, 0.0f), 1.0f) * 255.0f + 0.5f);
    return ru | (gu << 8) | (bu << 16);
}

static __forceinline__ __device__ int ClampInt(int v, int lo, int hi)
{
    return v < lo ? lo : (v > hi ? hi : v);
}

static __forceinline__ __device__ void TryWriteTextureFeedback(
    const LaunchParams::InstanceGeomRef &geomRef,
    unsigned int primitiveIndex,
    float u,
    float v,
    float uu,
    float vv,
    int textureWidth,
    int textureHeight,
    unsigned int mip);

static __forceinline__ __device__ float Clamp01(float v)
{
    return fminf(fmaxf(v, 0.0f), 1.0f);
}

static __forceinline__ __device__ float ApplyWrapMode(float uv, int wrapMode, bool *outBlack)
{
    *outBlack = false;
    switch (wrapMode)
    {
        case kWrapModeClamp:
            return Clamp01(uv);
        case kWrapModeMirror:
        {
            float t = fmodf(uv, 2.0f);
            if (t < 0.0f)
            {
                t += 2.0f;
            }
            const float mirrored = (t <= 1.0f) ? t : (2.0f - t);
            return Clamp01(mirrored);
        }
        case kWrapModeBlack:
            if (uv < 0.0f || uv > 1.0f)
            {
                *outBlack = true;
                return 0.0f;
            }
            return uv;
        case kWrapModeRepeat:
        case kWrapModeUseMetadata:
        case kWrapModeUnknown:
        default:
            return uv - floorf(uv);
    }
}

static __forceinline__ __device__ float3 MaterialSampleToViewColor(const float4 &sample)
{
    if (params.textureViewSemantic == kSemanticRoughness || params.textureViewSemantic == kSemanticMetallic ||
        params.textureViewSemantic == kSemanticOcclusion || params.textureViewSemantic == kSemanticIor ||
        params.textureViewSemantic == kSemanticOpacity)
    {
        return make_float3(sample.x, sample.x, sample.x);
    }
    return make_float3(sample.x, sample.y, sample.z);
}

static __forceinline__ __device__ bool TrySampleVirtualTexture(const LaunchParams::InstanceGeomRef &geomRef,
                                                               unsigned int primitiveIndex,
                                                               const LaunchParams::MaterialTextureRef &textureRef,
                                                               float u,
                                                               float v,
                                                               float uu,
                                                               float vv,
                                                               float3 &outColor)
{
    if (params.virtualTextureEnabled == 0)
    {
        return false;
    }

    if (params.virtualTextureTileEntries == 0ull || params.virtualTextureTilePixels == 0ull ||
        params.virtualTextureTileEntryCapacity <= 0)
    {
        outColor = make_float3(1.0f, 0.0f, 1.0f);
        return true;
    }

    bool blackS = false;
    bool blackT = false;
    const float wrappedU = ApplyWrapMode(uu, textureRef.wrapS, &blackS);
    const float wrappedV = ApplyWrapMode(vv, textureRef.wrapT, &blackT);
    if (blackS || blackT)
    {
        outColor = make_float3(0.0f, 0.0f, 0.0f);
        return true;
    }

    const int tileSize = max(params.feedbackTileSize, 1);
    const int texelX = ybi::texture::TexelFromUnitUV(wrappedU, textureRef.width);
    const int texelY = ybi::texture::TexelFromUnitUV(wrappedV, textureRef.height);
    const unsigned int tileX = ybi::texture::TileCoordFromTexel(texelX, tileSize);
    const unsigned int tileY = ybi::texture::TileCoordFromTexel(texelY, tileSize);
    const int udim = ybi::texture::UdimFromUV(u, v);
    const unsigned int udimBits = ybi::texture::UdimBitsFromUdim(udim);
    const unsigned int textureId =
        (unsigned int)ClampInt(geomRef.materialIndex, 0, int((1u << 23u) - 1u));
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
        outColor = make_float3(1.0f, 0.0f, 1.0f);
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
    const float4 sample = make_float4(float(pixels[sampleOffset + 0]) * (1.0f / 255.0f),
                                      float(pixels[sampleOffset + 1]) * (1.0f / 255.0f),
                                      float(pixels[sampleOffset + 2]) * (1.0f / 255.0f),
                                      float(pixels[sampleOffset + 3]) * (1.0f / 255.0f));
    outColor = MaterialSampleToViewColor(sample);
    TryWriteTextureFeedback(
        geomRef, primitiveIndex, u, v, wrappedU, wrappedV, textureRef.width, textureRef.height, mip);
    return true;
}

static __forceinline__ __device__ void TryWriteTextureFeedback(const LaunchParams::InstanceGeomRef &geomRef,
                                                               unsigned int primitiveIndex,
                                                               float u,
                                                               float v,
                                                               float uu,
                                                               float vv,
                                                               int textureWidth,
                                                               int textureHeight,
                                                               unsigned int mip)
{
    if (params.feedbackKeys == 0ull || params.feedbackStats == 0ull || params.feedbackCapacity <= 0 ||
        params.feedbackSamplePercent <= 0)
    {
        return;
    }

    const uint3 launchIndex = optixGetLaunchIndex();
    unsigned int seed =
        Hash32((launchIndex.x + 1u) * 73856093u ^ (launchIndex.y + 1u) * 19349663u ^
               (primitiveIndex + 1u) * 83492791u ^
               (static_cast<unsigned int>(params.currentSpp) + 1u) * 2654435761u ^
               (static_cast<unsigned int>(max(geomRef.materialIndex, 0)) + 1u) * 374761393u);
    if ((seed % 100u) >= (unsigned int)params.feedbackSamplePercent)
    {
        return;
    }

    const int tileSize = max(params.feedbackTileSize, 1);
    const int texelX = ybi::texture::TexelFromUnitUV(uu, textureWidth);
    const int texelY = ybi::texture::TexelFromUnitUV(vv, textureHeight);
    const unsigned int tileX = ybi::texture::TileCoordFromTexel(texelX, tileSize);
    const unsigned int tileY = ybi::texture::TileCoordFromTexel(texelY, tileSize);

    const int udim = ybi::texture::UdimFromUV(u, v);
    const unsigned int udimBits = ybi::texture::UdimBitsFromUdim(udim);

    const unsigned int textureId =
        (unsigned int)ClampInt(geomRef.materialIndex, 0, int((1u << 23u) - 1u));
    const unsigned long long key =
        ybi::texture::PackVirtualTextureKey(tileX, tileY, udimBits, textureId, mip);

    unsigned int *stats = reinterpret_cast<unsigned int *>(params.feedbackStats);
    unsigned long long *keys = reinterpret_cast<unsigned long long *>(params.feedbackKeys);
    const unsigned int index = atomicAdd(&stats[0], 1u);
    if (index < (unsigned int)params.feedbackCapacity)
    {
        keys[index] = key;
    }
    else
    {
        atomicAdd(&stats[1], 1u);
    }
}

static __forceinline__ __device__ bool
TryWriteFeedbackOnly(const LaunchParams::InstanceGeomRef &geomRef,
                     unsigned int primitiveIndex,
                     const float3 &barycentrics,
                     unsigned int mip)
{
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

    const float2_simple *texcoords = reinterpret_cast<const float2_simple *>(geomRef.texcoords);
    const float2_simple uv0 = texcoords[t0];
    const float2_simple uv1 = texcoords[t1];
    const float2_simple uv2 = texcoords[t2];
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
    if (params.textureViewSemantic < 0 || params.textureViewSemantic >= params.materialTextureRefSemanticCount)
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
        geomRef, primitiveIndex, u, v, uu, vv, textureRef.width, textureRef.height, mip);
    return true;
}

static __forceinline__ __device__ float3 SkyColor(const float3 &direction)
{
    const float t = 0.5f * (direction.y + 1.0f);
    const float3 top = make_float3(0.7f, 0.8f, 1.0f);
    const float3 bottom = make_float3(0.2f, 0.25f, 0.35f);
    return make_float3((1.0f - t) * top.x + t * bottom.x,
                       (1.0f - t) * top.y + t * bottom.y,
                       (1.0f - t) * top.z + t * bottom.z);
}

static __forceinline__ __device__ bool
TrySampleMaterialTexture(const LaunchParams::InstanceGeomRef &geomRef,
                         unsigned int primitiveIndex,
                         const float3 &barycentrics,
                         float3 &outColor)
{
    if (params.materialTextureRefs == 0ull || params.materialTextureRefCount <= 0 ||
        params.materialTextureRefStride <= 0)
    {
        const unsigned int prev = atomicOr(&g_try_sample_fail_mask, 1u << 0);
        if ((prev & (1u << 0)) == 0u)
        {
            printf("TrySampleMaterialTexture fail[0]: material ref invalid ptr=%llu count=%d stride=%d mat=%d\n",
                   params.materialTextureRefs,
                   params.materialTextureRefCount,
                   params.materialTextureRefStride,
                   geomRef.materialIndex);
        }
        return false;
    }
    if (geomRef.texcoords == 0ull || geomRef.texcoordIndices == 0ull)
    {
        const unsigned int prev = atomicOr(&g_try_sample_fail_mask, 1u << 1);
        if ((prev & (1u << 1)) == 0u)
        {
            printf("TrySampleMaterialTexture fail[1]: texcoord buffers missing tc=%llu tci=%llu nTc=%d nTci=%d\n",
                   geomRef.texcoords,
                   geomRef.texcoordIndices,
                   geomRef.numTexcoords,
                   geomRef.numTexcoordIndices);
        }
        return false;
    }

    const int triCornerBase = int(primitiveIndex) * 3;
    if (triCornerBase + 2 >= geomRef.numTexcoordIndices)
    {
        const unsigned int prev = atomicOr(&g_try_sample_fail_mask, 1u << 2);
        if ((prev & (1u << 2)) == 0u)
        {
            printf("TrySampleMaterialTexture fail[2]: tri texcoord range prim=%u base=%d count=%d\n",
                   primitiveIndex,
                   triCornerBase,
                   geomRef.numTexcoordIndices);
        }
        return false;
    }

    const int *tcIndices = reinterpret_cast<const int *>(geomRef.texcoordIndices);
    const int t0 = tcIndices[triCornerBase + 0];
    const int t1 = tcIndices[triCornerBase + 1];
    const int t2 = tcIndices[triCornerBase + 2];
    if (t0 < 0 || t0 >= geomRef.numTexcoords || t1 < 0 || t1 >= geomRef.numTexcoords || t2 < 0 ||
        t2 >= geomRef.numTexcoords)
    {
        const unsigned int prev = atomicOr(&g_try_sample_fail_mask, 1u << 3);
        if ((prev & (1u << 3)) == 0u)
        {
            printf("TrySampleMaterialTexture fail[3]: tc index out of range (%d,%d,%d) nTc=%d\n",
                   t0,
                   t1,
                   t2,
                   geomRef.numTexcoords);
        }
        return false;
    }

    const float2_simple *texcoords = reinterpret_cast<const float2_simple *>(geomRef.texcoords);
    const float2_simple uv0 = texcoords[t0];
    const float2_simple uv1 = texcoords[t1];
    const float2_simple uv2 = texcoords[t2];
    const float u = uv0.x * barycentrics.x + uv1.x * barycentrics.y + uv2.x * barycentrics.z;
    const float v = uv0.y * barycentrics.x + uv1.y * barycentrics.y + uv2.y * barycentrics.z;
    const int udim = ybi::texture::UdimFromUV(u, v);
    const int udimSlot = ybi::texture::UdimSlotFromUdim(udim, params.materialTextureRefStride);
    const float uu = ybi::texture::UnitUV(u);
    const float vv = ybi::texture::UnitUV(v);

    if (geomRef.materialIndex < 0 || geomRef.materialIndex >= params.materialTextureRefCount)
    {
        const unsigned int prev = atomicOr(&g_try_sample_fail_mask, 1u << 5);
        if ((prev & (1u << 5)) == 0u)
        {
            printf("TrySampleMaterialTexture fail[5]: material index invalid mat=%d count=%d\n",
                   geomRef.materialIndex,
                   params.materialTextureRefCount);
        }
        return false;
    }
    const int materialIndex = geomRef.materialIndex;
    if (params.textureViewSemantic < 0 || params.textureViewSemantic >= params.materialTextureRefSemanticCount)
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
        const unsigned int prev = atomicOr(&g_try_sample_fail_mask, 1u << 4);
        if ((prev & (1u << 4)) == 0u)
        {
            printf("TrySampleMaterialTexture fail[4]: udim slot invalid slot=%d max=%d mat=%d udim=%d\n",
                   slot,
                   maxSlots,
                   materialIndex,
                   udim);
        }
        return false;
    }

    LaunchParams::MaterialTextureRef textureRef = materialRefs[slot];
    if (params.virtualTextureEnabled != 0)
    {
        if (TrySampleVirtualTexture(
                geomRef, primitiveIndex, textureRef, u, v, uu, vv, outColor))
        {
            return true;
        }
    }

    if (textureRef.textureObject == 0ull || textureRef.valid == 0)
    {
        const unsigned int prev = atomicOr(&g_try_sample_fail_mask, 1u << 4);
        if ((prev & (1u << 4)) == 0u)
        {
            printf("TrySampleMaterialTexture fail[4]: texture ref invalid obj=%llu valid=%d mat=%d udim=%d\n",
                   textureRef.textureObject,
                   textureRef.valid,
                   materialIndex,
                   udim);
        }
        return false;
    }

    const cudaTextureObject_t textureObject =
        static_cast<cudaTextureObject_t>(textureRef.textureObject);
    const float4 sample = tex2D<float4>(textureObject, uu, vv);
    outColor = MaterialSampleToViewColor(sample);
    TryWriteTextureFeedback(
        geomRef, primitiveIndex, u, v, uu, vv, textureRef.width, textureRef.height, 0u);
    return true;
}

static __forceinline__ __device__ float3 ComputeDirection(const uint3 &launchIndex,
                                                          const uint3 &launchDims,
                                                          const float2 &pixelOffset)
{
    const float2 ndc =
        make_float2(((float)launchIndex.x + pixelOffset.x) / (float)launchDims.x * 2.0f - 1.0f,
                    1.0f - ((float)launchIndex.y + pixelOffset.y) / (float)launchDims.y * 2.0f);

    return Normalize3(
        make_float3(params.cameraU.x * ndc.x + params.cameraV.x * ndc.y + params.cameraW.x,
                    params.cameraU.y * ndc.x + params.cameraV.y * ndc.y + params.cameraW.y,
                    params.cameraU.z * ndc.x + params.cameraV.z * ndc.y + params.cameraW.z));
}

static __forceinline__ __device__ unsigned int TraceColor(const float3 &origin,
                                                          const float3 &direction)
{
    unsigned int packedColor = 0;
    optixTrace(params.traversable,
               origin,
               direction,
               0.001f,
               1e20f,
               0.0f,
               OptixVisibilityMask(255),
               OPTIX_RAY_FLAG_NONE,
               0,
               1,
               0,
               packedColor);
    return packedColor;
}

static __forceinline__ __device__ unsigned int
TraceOcclusion(const float3 &origin, const float3 &direction, float tMin, float tMax)
{
    unsigned int hit = 0x80000000u;
    optixTrace(params.traversable,
               origin,
               direction,
               tMin,
               tMax,
               0.0f,
               OptixVisibilityMask(255),
               OPTIX_RAY_FLAG_DISABLE_CLOSESTHIT | OPTIX_RAY_FLAG_ENFORCE_ANYHIT |
                   OPTIX_RAY_FLAG_TERMINATE_ON_FIRST_HIT,
               0,
               1,
               0,
               hit);
    return hit;
}

extern "C" __global__ void __raygen__primary()
{
    const uint3 launchIndex = optixGetLaunchIndex();
    const uint3 launchDims = optixGetLaunchDimensions();
    const unsigned int pixelIndex = launchIndex.y * launchDims.x + launchIndex.x;

    const float2 centerOffset = make_float2(0.5f, 0.5f);
    const float3 centerDirection = ComputeDirection(launchIndex, launchDims, centerOffset);
    const unsigned int packedColor = TraceColor(params.cameraOrigin, centerDirection);
    uchar4 *image = reinterpret_cast<uchar4 *>(params.image);
    image[pixelIndex] = make_uchar4((unsigned char)(packedColor & 255u),
                                    (unsigned char)((packedColor >> 8) & 255u),
                                    (unsigned char)((packedColor >> 16) & 255u),
                                    255u);
}

extern "C" __global__ void __raygen__feedback()
{
    const uint3 launchIndex = optixGetLaunchIndex();
    const uint3 launchDims = optixGetLaunchDimensions();
    const float2 centerOffset = make_float2(0.5f, 0.5f);
    const float3 centerDirection = ComputeDirection(launchIndex, launchDims, centerOffset);
    TraceColor(params.cameraOrigin, centerDirection);
}

extern "C" __global__ void __miss__primary()
{
    const unsigned int payload = optixGetPayload_0();
    if (payload == 0x80000000u)
    {
        optixSetPayload_0(0u);
        return;
    }

    const float3 direction = optixGetWorldRayDirection();
    const float3 color = SkyColor(direction);
    optixSetPayload_0(PackColor(color.x, color.y, color.z));
}

extern "C" __global__ void __anyhit__primary()
{
    const unsigned int payload = optixGetPayload_0();
    if (payload == 0x80000000u)
    {
        optixSetPayload_0(1u);
        optixTerminateRay();
    }
}

extern "C" __global__ void __closesthit__primary()
{
    if (params.integrator == 2)
    {
        const unsigned int hitKind = optixGetHitKind();
        if (hitKind == OPTIX_HIT_KIND_TRIANGLE_FRONT_FACE ||
            hitKind == OPTIX_HIT_KIND_TRIANGLE_BACK_FACE)
        {
            const float2 bary = optixGetTriangleBarycentrics();
            const float3 barycentrics = make_float3(1.0f - bary.x - bary.y, bary.x, bary.y);
            const unsigned int instanceId = optixGetInstanceId();
            if (params.instanceGeomRefs != 0ull &&
                instanceId < (unsigned int)params.instanceGeomRefCount)
            {
                const LaunchParams::InstanceGeomRef *refs =
                    reinterpret_cast<const LaunchParams::InstanceGeomRef *>(params.instanceGeomRefs);
                const LaunchParams::InstanceGeomRef ref = refs[instanceId];
                TryWriteFeedbackOnly(ref, optixGetPrimitiveIndex(), barycentrics, 1u);
            }
        }
        optixSetPayload_0(0u);
        return;
    }

    if (params.integrator == 1)
    {
        const int spp = max(params.spp, 1);
        const uint3 launchIndex = optixGetLaunchIndex();
        unsigned int rngState = Hash32((launchIndex.x + 1u) * 73856093u ^
                                       (launchIndex.y + 1u) * 19349663u ^
                                       (static_cast<unsigned int>(params.currentSpp) + 1u) * 83492791u);

        const float3 rayOrigin = optixGetWorldRayOrigin();
        const float3 rayDirection = Normalize3(optixGetWorldRayDirection());
        const float tHit = optixGetRayTmax();
        const float3 hitPoint = make_float3(rayOrigin.x + rayDirection.x * tHit,
                                            rayOrigin.y + rayDirection.y * tHit,
                                            rayOrigin.z + rayDirection.z * tHit);
        float3 normal = Normalize3(make_float3(-rayDirection.x, -rayDirection.y, -rayDirection.z));
        const unsigned int hitKind = optixGetHitKind();
        if (hitKind == OPTIX_HIT_KIND_TRIANGLE_FRONT_FACE ||
            hitKind == OPTIX_HIT_KIND_TRIANGLE_BACK_FACE)
        {
            const unsigned int instanceId = optixGetInstanceId();
            if (params.instanceGeomRefs != 0ull &&
                instanceId < (unsigned int)params.instanceGeomRefCount)
            {
                const LaunchParams::InstanceGeomRef *refs =
                    reinterpret_cast<const LaunchParams::InstanceGeomRef *>(
                        params.instanceGeomRefs);
                const LaunchParams::InstanceGeomRef ref = refs[instanceId];
                const int primitiveIndex = int(optixGetPrimitiveIndex());
                const int indexBase = primitiveIndex * 3;
                const float3 *positions = reinterpret_cast<const float3 *>(ref.positions);
                const int *indices = reinterpret_cast<const int *>(ref.indices);
                if (positions != nullptr && indices != nullptr && indexBase + 2 < ref.numIndices)
                {
                    const int i0 = indices[indexBase + 0];
                    const int i1 = indices[indexBase + 1];
                    const int i2 = indices[indexBase + 2];
                    if (i0 >= 0 && i0 < ref.numPositions && i1 >= 0 && i1 < ref.numPositions &&
                        i2 >= 0 && i2 < ref.numPositions)
                    {
                        const float3 p0 = positions[i0];
                        const float3 p1 = positions[i1];
                        const float3 p2 = positions[i2];
                        const float3 edge01 = make_float3(p1.x - p0.x, p1.y - p0.y, p1.z - p0.z);
                        const float3 edge02 = make_float3(p2.x - p0.x, p2.y - p0.y, p2.z - p0.z);
                        const float3 localNormal = Normalize3(Cross3(edge01, edge02));
                        const float3 worldNormal =
                            Normalize3(optixTransformNormalFromObjectToWorldSpace(localNormal));
                        normal = FaceForward(worldNormal, rayDirection);
                    }
                }
            }
        }
        float3 tangent, bitangent;
        BuildOrthonormalBasis(normal, tangent, bitangent);

        int visible = 0;
        for (int i = 0; i < spp; i++)
        {
            const float u1 = Random01(rngState);
            const float u2 = Random01(rngState);
            const float3 local = SampleCosineHemisphere(u1, u2);
            const float3 sampleDir = Normalize3(
                make_float3(tangent.x * local.x + bitangent.x * local.y + normal.x * local.z,
                            tangent.y * local.x + bitangent.y * local.y + normal.y * local.z,
                            tangent.z * local.x + bitangent.z * local.y + normal.z * local.z));
            const float3 sampleOrigin = make_float3(hitPoint.x + normal.x * params.aoBias,
                                                    hitPoint.y + normal.y * params.aoBias,
                                                    hitPoint.z + normal.z * params.aoBias);
            const unsigned int occluded =
                TraceOcclusion(sampleOrigin, sampleDir, params.aoBias, params.aoMaxDistance);
            if (occluded == 0u)
            {
                visible++;
            }
        }

        const float ao = float(visible) / float(spp);
        const float3 aoColor = make_float3(ao, ao, ao);
        optixSetPayload_0(PackColor(aoColor.x, aoColor.y, aoColor.z));
        return;
    }

    float3 color = make_float3(0.7f, 0.7f, 0.7f);
    const unsigned int hitKind = optixGetHitKind();
    if (hitKind == OPTIX_HIT_KIND_TRIANGLE_FRONT_FACE ||
        hitKind == OPTIX_HIT_KIND_TRIANGLE_BACK_FACE)
    {
        const float2 bary = optixGetTriangleBarycentrics();
        const float3 barycentrics = make_float3(1.0f - bary.x - bary.y, bary.x, bary.y);
        const unsigned int instanceId = optixGetInstanceId();
        if (params.instanceGeomRefs != 0ull &&
            instanceId < (unsigned int)params.instanceGeomRefCount)
        {
            const LaunchParams::InstanceGeomRef *refs =
                reinterpret_cast<const LaunchParams::InstanceGeomRef *>(params.instanceGeomRefs);
            const LaunchParams::InstanceGeomRef ref = refs[instanceId];
            const bool sampled =
                TrySampleMaterialTexture(ref, optixGetPrimitiveIndex(), barycentrics, color);
            if (!sampled)
            {
                color = make_float3(0.0f, 0.0f, 0.0f);
            }
            else if (atomicCAS(&g_try_sample_logged, 0u, 1u) == 0u)
            {
                printf("TrySampleMaterialTexture succeeded\n");
            }
        }
    }
    else
    {
        color = make_float3(0.7f, 0.7f, 0.7f);
    }
    optixSetPayload_0(PackColor(color.x, color.y, color.z));
}
