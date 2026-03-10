#pragma once

#ifndef YBI_NAMESPACE_BEGIN
#define YBI_NAMESPACE_BEGIN namespace ybi {
#endif

#ifndef YBI_NAMESPACE_END
#define YBI_NAMESPACE_END }
#endif

#include "device/device.h"
#include "scene/material_light.h"
#include "util/vec3.h"

YBI_NAMESPACE_BEGIN

struct LaunchParams
{
    struct InstanceGeomRef
    {
        unsigned long long positions;
        unsigned long long indices;
        unsigned long long texcoords;
        unsigned long long texcoordIndices;
        unsigned long long normals;
        unsigned long long normalIndices;
        int numPositions;
        int numIndices;
        int numTexcoords;
        int numTexcoordIndices;
        int numNormals;
        int numNormalIndices;
        int materialIndex;
        int _padding0;
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

    struct VirtualTextureMipInfo
    {
        unsigned int level;
        unsigned int udimInfoOffset;
        unsigned int udimInfoCount;
        unsigned int _padding0;
    };

    struct VirtualTextureUdimInfo
    {
        unsigned int basePageX;
        unsigned int basePageY;
        unsigned int pageCountX;
        unsigned int pageCountY;
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
    ybi::Vec3 cameraOrigin;
    ybi::Vec3 cameraU;
    ybi::Vec3 cameraV;
    ybi::Vec3 cameraW;
    WireframeConfig wireframe;
    int integrator;
    int spp;
    int maxDepth;
    float aoBias;
    float aoMaxDistance;
    unsigned long long instanceGeomRefs;
    int instanceGeomRefCount;
    unsigned long long materialTextureRefs;
    int materialTextureRefCount;
    int materialTextureRefStride;
    int materialTextureRefSemanticCount;
    MaterialTextureRef domeTextureRef;
    unsigned long long materials;
    int materialCount;
    unsigned long long lights;
    int lightCount;
    unsigned long long lightShadowExcludeRefs;
    int lightShadowExcludeRefCount;
    int textureViewSemantic;
    unsigned long long feedbackKeys;
    unsigned long long feedbackStats;
    int feedbackCapacity;
    int feedbackSamplePercent;
    int feedbackTileSize;
    int currentSpp;
    int singlePixelEnabled;
    int singlePixelX;
    int singlePixelY;
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
    unsigned long long virtualTextureUdimInfos;
    int virtualTextureUdimInfoCount;
};

YBI_NAMESPACE_END
