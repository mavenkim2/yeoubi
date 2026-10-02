#pragma once

#include "core/affine_transform.h"
#include "core/vec2f.h"
#include "core/vec3f.h"
#include "device/cuda_compat.h"
#include "kernel/integrator/rng.h"

namespace ybi
{

enum KernelType
{
    KERNEL_PRIMARY_RAY,
    KERNEL_SHADE_HIT,
};

enum IntersectionType
{
    INTERSECTION_TYPE_NONE = 0,
    INTERSECTION_TYPE_TRIANGLE = 1,
};

struct IntersectionInfo
{
    float t;
    Vec2f bary;
    AffineTransform objectToWorld;
    uint32_t primIndex;
    uint32_t clusterId;
    uint32_t id;
    IntersectionType type;
};

struct RayCone
{
    float width;
    float spreadAngle;
};

struct RNG;

struct WavefrontSampleState
{
    uint32_t *pixelX;
    uint32_t *pixelY;
    RNG *rng;
    RayCone *rayCone;
    Vec3f *radiance;
    Vec3f *throughput;
    Vec3f *previousSurfacePosition;
    float *previousBsdfPdf;
    Vec3f *rayOrigin;
    Vec3f *rayDirection;
    uint32_t *depth;
    KernelType *kernelTypes;

    // Intersection state
    float *t;
    Vec2f *bary;
    AffineTransform *objectToWorld;
    uint32_t *primIndex;
    uint32_t *clusterId;
    uint32_t *objectId;
    IntersectionType *type;

    YBI_GPU uint32_t ReadPixelX(uint32_t sampleIndex) const
    {
        return pixelX[sampleIndex];
    }

    YBI_GPU void WritePixelX(uint32_t sampleIndex, uint32_t value) const
    {
        pixelX[sampleIndex] = value;
    }

    YBI_GPU uint32_t ReadPixelY(uint32_t sampleIndex) const
    {
        return pixelY[sampleIndex];
    }

    YBI_GPU void WritePixelY(uint32_t sampleIndex, uint32_t value) const
    {
        pixelY[sampleIndex] = value;
    }

    YBI_GPU RNG ReadRNG(uint32_t sampleIndex) const
    {
        return rng[sampleIndex];
    }

    YBI_GPU void WriteRNG(uint32_t sampleIndex, const RNG &value) const
    {
        rng[sampleIndex] = value;
    }

    YBI_GPU RayCone ReadRayCone(uint32_t sampleIndex) const
    {
        return rayCone[sampleIndex];
    }

    YBI_GPU void WriteRayCone(uint32_t sampleIndex, const RayCone &value) const
    {
        rayCone[sampleIndex] = value;
    }

    YBI_GPU Vec3f ReadRadiance(uint32_t sampleIndex) const
    {
        return radiance[sampleIndex];
    }

    YBI_GPU void WriteRadiance(uint32_t sampleIndex, const Vec3f &value) const
    {
        radiance[sampleIndex] = value;
    }

    YBI_GPU Vec3f ReadThroughput(uint32_t sampleIndex) const
    {
        return throughput[sampleIndex];
    }

    YBI_GPU void WriteThroughput(uint32_t sampleIndex, const Vec3f &value) const
    {
        throughput[sampleIndex] = value;
    }

    YBI_GPU Vec3f ReadPreviousSurfacePosition(uint32_t sampleIndex) const
    {
        return previousSurfacePosition[sampleIndex];
    }

    YBI_GPU void WritePreviousSurfacePosition(uint32_t sampleIndex, const Vec3f &value) const
    {
        previousSurfacePosition[sampleIndex] = value;
    }

    YBI_GPU float ReadPreviousBsdfPdf(uint32_t sampleIndex) const
    {
        return previousBsdfPdf[sampleIndex];
    }

    YBI_GPU void WritePreviousBsdfPdf(uint32_t sampleIndex, float value) const
    {
        previousBsdfPdf[sampleIndex] = value;
    }

    YBI_GPU Vec3f ReadRayOrigin(uint32_t sampleIndex) const
    {
        return rayOrigin[sampleIndex];
    }

    YBI_GPU void WriteRayOrigin(uint32_t sampleIndex, const Vec3f &value) const
    {
        rayOrigin[sampleIndex] = value;
    }

    YBI_GPU Vec3f ReadRayDirection(uint32_t sampleIndex) const
    {
        return rayDirection[sampleIndex];
    }

    YBI_GPU void WriteRayDirection(uint32_t sampleIndex, const Vec3f &value) const
    {
        rayDirection[sampleIndex] = value;
    }

    YBI_GPU uint32_t ReadDepth(uint32_t sampleIndex) const
    {
        return depth[sampleIndex];
    }

    YBI_GPU void WriteDepth(uint32_t sampleIndex, uint32_t value) const
    {
        depth[sampleIndex] = value;
    }

    YBI_GPU KernelType ReadKernelType(uint32_t sampleIndex) const
    {
        return kernelTypes[sampleIndex];
    }

    YBI_GPU void WriteKernelType(uint32_t sampleIndex, KernelType value) const
    {
        kernelTypes[sampleIndex] = value;
    }

    YBI_GPU float ReadT(uint32_t sampleIndex) const
    {
        return t[sampleIndex];
    }

    YBI_GPU void WriteT(uint32_t sampleIndex, float value) const
    {
        t[sampleIndex] = value;
    }

    YBI_GPU Vec2f ReadBary(uint32_t sampleIndex) const
    {
        return bary[sampleIndex];
    }

    YBI_GPU void WriteBary(uint32_t sampleIndex, const Vec2f &value) const
    {
        bary[sampleIndex] = value;
    }

    YBI_GPU AffineTransform ReadObjectToWorld(uint32_t sampleIndex) const
    {
        return objectToWorld[sampleIndex];
    }

    YBI_GPU void WriteObjectToWorld(uint32_t sampleIndex, const AffineTransform &value) const
    {
        objectToWorld[sampleIndex] = value;
    }

    YBI_GPU uint32_t ReadPrimIndex(uint32_t sampleIndex) const
    {
        return primIndex[sampleIndex];
    }

    YBI_GPU void WritePrimIndex(uint32_t sampleIndex, uint32_t value) const
    {
        primIndex[sampleIndex] = value;
    }

    YBI_GPU uint32_t ReadClusterId(uint32_t sampleIndex) const
    {
        return clusterId[sampleIndex];
    }

    YBI_GPU void WriteClusterId(uint32_t sampleIndex, uint32_t value) const
    {
        clusterId[sampleIndex] = value;
    }

    YBI_GPU uint32_t ReadObjectId(uint32_t sampleIndex) const
    {
        return objectId[sampleIndex];
    }

    YBI_GPU void WriteId(uint32_t sampleIndex, uint32_t value) const
    {
        objectId[sampleIndex] = value;
    }

    YBI_GPU IntersectionType ReadType(uint32_t sampleIndex) const
    {
        return type[sampleIndex];
    }

    YBI_GPU void WriteType(uint32_t sampleIndex, IntersectionType value) const
    {
        type[sampleIndex] = value;
    }

    YBI_GPU IntersectionInfo ReadIntersection(uint32_t sampleIndex) const
    {
        IntersectionInfo info;
        info.t = ReadT(sampleIndex);
        info.bary = ReadBary(sampleIndex);
        info.objectToWorld = ReadObjectToWorld(sampleIndex);
        info.primIndex = ReadPrimIndex(sampleIndex);
        info.clusterId = ReadClusterId(sampleIndex);
        info.id = ReadObjectId(sampleIndex);
        info.type = ReadType(sampleIndex);
        return info;
    }

    YBI_GPU void WriteIntersection(uint32_t sampleIndex, const IntersectionInfo &info) const
    {
        WriteT(sampleIndex, info.t);
        WriteBary(sampleIndex, info.bary);
        WriteObjectToWorld(sampleIndex, info.objectToWorld);
        WritePrimIndex(sampleIndex, info.primIndex);
        WriteClusterId(sampleIndex, info.clusterId);
        WriteId(sampleIndex, info.id);
        WriteType(sampleIndex, info.type);
    }
};

struct SampleState
{
    uint32_t pixelX;
    uint32_t pixelY;
    RNG rng;
    RayCone rayCone;

    Vec3f radiance;
    Vec3f throughput;
    Vec3f previousSurfacePosition;
    float previousBsdfPdf;

    Vec3f rayOrigin;
    Vec3f rayDirection;

    uint32_t depth;
};

} // namespace ybi
