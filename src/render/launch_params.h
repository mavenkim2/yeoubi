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
};

YBI_NAMESPACE_END
