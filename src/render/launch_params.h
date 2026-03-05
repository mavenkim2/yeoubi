#pragma once

#ifndef YBI_NAMESPACE_BEGIN
#define YBI_NAMESPACE_BEGIN namespace ybi {
#endif

#ifndef YBI_NAMESPACE_END
#define YBI_NAMESPACE_END }
#endif

#include "device/device.h"
#include "util/float3.h"

YBI_NAMESPACE_BEGIN

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

    struct VirtualTextureMipInfo
    {
        unsigned int level;
        unsigned int basePageX;
        unsigned int basePageY;
        unsigned int pageCountX;
        unsigned int pageCountY;
        unsigned int pagesPerUdimX;
        unsigned int pagesPerUdimY;
        unsigned int _padding0;
    };

    struct VirtualTextureTextureMeta
    {
        unsigned int mipInfoOffset;
        unsigned int mipCount;
        unsigned int tailFirstMip;
        unsigned int activeUdimCount;
        unsigned long long tailPixels;
        unsigned int tailPageCountX;
        unsigned int tailPageCountY;
        unsigned int _padding0;
        short udimToLocal[128];
    };

    BVHHandle traversable;
    DevicePtr image;
    int width;
    int height;
    ybi::float3 cameraOrigin;
    ybi::float3 cameraU;
    ybi::float3 cameraV;
    ybi::float3 cameraW;
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

    unsigned long long virtualTexturePageTableEntries;
    unsigned long long virtualTexturePageTableMipOffsets;
    unsigned long long virtualTexturePageTableMipWidths;
    unsigned long long virtualTexturePageTableMipHeights;
    int virtualTexturePageTableMipCount;
    int virtualTexturePageSize;
    unsigned long long virtualTextureStreamPixels;
    int virtualTextureStreamPageCountX;
    int virtualTextureStreamPageCountY;
    int virtualTextureSampleMip;
    unsigned long long virtualTextureTextureMeta;
    int virtualTextureTextureMetaCount;
    unsigned long long virtualTextureMipInfos;
    int virtualTextureMipInfoCount;
};

YBI_NAMESPACE_END
