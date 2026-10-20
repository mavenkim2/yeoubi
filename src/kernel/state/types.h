#pragma once

#include "core/affine_transform.h"
#include "core/vec2f.h"
#include "core/vec3f.h"
#include "device/cuda_compat.h"
#include "kernel/bsdf/closure_types.h"
#include "kernel/sampling/rng.h"

namespace ybi
{

enum KernelType
{
    KERNEL_GENERATE_CAMERA_RAYS,
    KERNEL_INTERSECT,
    KERNEL_SHADE_HIT,
    KERNEL_DIRECT_LIGHTING,
    KERNEL_INDIRECT_LIGHTING,
    KERNEL_FINISHED,
    KERNEL_TYPE_COUNT,
};

struct KernelCounts
{
    uint32_t counts[KERNEL_TYPE_COUNT] = {};
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

struct TriangleSurfaceInfo
{
    Vec3f position;
    Vec3f geometricNormal;
    Vec3f shadingNormal;
    Vec3f dpdu;
    Vec2f uv;
    float triangleArea;
    float uvArea;
    uint32_t materialIndex;
};

struct RNG;

struct WavefrontSampleState
{
    uint32_t *pixelIndex;
    RNG *rng;
    RayCone *rayCone;
    Vec3f *radiance;
#if defined(YBI_USE_PATH_GUIDING)
    int *radianceRecordIndex;
#endif
    Vec3f *throughput;
    Vec3f *previousSurfacePosition;
    float *previousBsdfPdf;
    uint32_t *previousBsdfFlags;
    float *etaScale;
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

    // Surface state
    uint32_t *materialIndex;
    Vec3f *position;
    Vec3f *geometricNormal;
    Vec3f *shadingNormal;
    Vec3f *dpdu;
    Vec2f *uv;
    float *triangleArea;
    float *uvArea;
    float *closureData;

    YBI_GPU uint32_t ReadPixelIndex(uint32_t sampleIndex) const
    {
        return pixelIndex[sampleIndex];
    }

    YBI_GPU void WritePixelIndex(uint32_t sampleIndex, uint32_t value) const
    {
        pixelIndex[sampleIndex] = value;
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

#if defined(YBI_USE_PATH_GUIDING)
    YBI_GPU int ReadRadianceRecordIndex(uint32_t sampleIndex) const
    {
        return radianceRecordIndex[sampleIndex];
    }

    YBI_GPU void WriteRadianceRecordIndex(uint32_t sampleIndex, int value) const
    {
        radianceRecordIndex[sampleIndex] = value;
    }
#endif

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

    YBI_GPU uint32_t ReadPreviousBsdfFlags(uint32_t sampleIndex) const
    {
        return previousBsdfFlags[sampleIndex];
    }

    YBI_GPU void WritePreviousBsdfFlags(uint32_t sampleIndex, uint32_t value) const
    {
        previousBsdfFlags[sampleIndex] = value;
    }

    YBI_GPU float ReadEtaScale(uint32_t sampleIndex) const
    {
        return etaScale[sampleIndex];
    }

    YBI_GPU void WriteEtaScale(uint32_t sampleIndex, float value) const
    {
        etaScale[sampleIndex] = value;
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

    YBI_GPU uint32_t ReadMaterialIndex(uint32_t sampleIndex) const
    {
        return materialIndex[sampleIndex];
    }

    YBI_GPU void WriteMaterialIndex(uint32_t sampleIndex, uint32_t value) const
    {
        materialIndex[sampleIndex] = value;
    }

    YBI_GPU Vec3f ReadPosition(uint32_t sampleIndex) const
    {
        return position[sampleIndex];
    }

    YBI_GPU void WritePosition(uint32_t sampleIndex, const Vec3f &value) const
    {
        position[sampleIndex] = value;
    }

    YBI_GPU Vec3f ReadGeometricNormal(uint32_t sampleIndex) const
    {
        return geometricNormal[sampleIndex];
    }

    YBI_GPU void WriteGeometricNormal(uint32_t sampleIndex, const Vec3f &value) const
    {
        geometricNormal[sampleIndex] = value;
    }

    YBI_GPU Vec3f ReadShadingNormal(uint32_t sampleIndex) const
    {
        return shadingNormal[sampleIndex];
    }

    YBI_GPU void WriteShadingNormal(uint32_t sampleIndex, const Vec3f &value) const
    {
        shadingNormal[sampleIndex] = value;
    }

    YBI_GPU Vec3f ReadDpdu(uint32_t sampleIndex) const
    {
        return dpdu[sampleIndex];
    }

    YBI_GPU void WriteDpdu(uint32_t sampleIndex, const Vec3f &value) const
    {
        dpdu[sampleIndex] = value;
    }

    YBI_GPU Vec2f ReadUV(uint32_t sampleIndex) const
    {
        return uv[sampleIndex];
    }

    YBI_GPU void WriteUV(uint32_t sampleIndex, const Vec2f &value) const
    {
        uv[sampleIndex] = value;
    }

    YBI_GPU float ReadTriangleArea(uint32_t sampleIndex) const
    {
        return triangleArea[sampleIndex];
    }

    YBI_GPU void WriteTriangleArea(uint32_t sampleIndex, float value) const
    {
        triangleArea[sampleIndex] = value;
    }

    YBI_GPU float ReadUVArea(uint32_t sampleIndex) const
    {
        return uvArea[sampleIndex];
    }

    YBI_GPU void WriteUVArea(uint32_t sampleIndex, float value) const
    {
        uvArea[sampleIndex] = value;
    }

    YBI_GPU ClosureInitializationInfo ReadClosureInitializationInfo(uint32_t sampleIndex) const
    {
        ClosureInitializationInfo info;
        info.geometricNormal = ReadGeometricNormal(sampleIndex);
        info.shadingNormal = ReadShadingNormal(sampleIndex);
        info.dpdu = ReadDpdu(sampleIndex);
        info.uv = ReadUV(sampleIndex);
        info.rayDirection = ReadRayDirection(sampleIndex);
        info.rayConeWidth = ReadRayCone(sampleIndex).width;
        info.triangleArea = ReadTriangleArea(sampleIndex);
        info.uvArea = ReadUVArea(sampleIndex);
        info.closureData = closureData + sampleIndex * MAX_CLOSURE_SIZE;
        return info;
    }

    YBI_GPU ClosureEvaluationInfo ReadClosureEvaluationInfo(uint32_t sampleIndex) const
    {
        ClosureEvaluationInfo info;
        info.geometricNormal = ReadGeometricNormal(sampleIndex);
        info.shadingNormal = ReadShadingNormal(sampleIndex);
        info.dpdu = ReadDpdu(sampleIndex);
        info.closureData = closureData + sampleIndex * MAX_CLOSURE_SIZE;
        return info;
    }

    YBI_GPU TriangleSurfaceInfo ReadTriangleInfo(uint32_t sampleIndex) const
    {
        TriangleSurfaceInfo info;
        info.materialIndex = ReadMaterialIndex(sampleIndex);
        info.position = ReadPosition(sampleIndex);
        info.geometricNormal = ReadGeometricNormal(sampleIndex);
        info.shadingNormal = ReadShadingNormal(sampleIndex);
        info.dpdu = ReadDpdu(sampleIndex);
        info.uv = ReadUV(sampleIndex);
        info.triangleArea = ReadTriangleArea(sampleIndex);
        info.uvArea = ReadUVArea(sampleIndex);
        return info;
    }

    YBI_GPU void WriteTriangleInfo(uint32_t sampleIndex, const TriangleSurfaceInfo &info) const
    {
        WriteMaterialIndex(sampleIndex, info.materialIndex);
        WritePosition(sampleIndex, info.position);
        WriteGeometricNormal(sampleIndex, info.geometricNormal);
        WriteShadingNormal(sampleIndex, info.shadingNormal);
        WriteDpdu(sampleIndex, info.dpdu);
        WriteUV(sampleIndex, info.uv);
        WriteTriangleArea(sampleIndex, info.triangleArea);
        WriteUVArea(sampleIndex, info.uvArea);
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
