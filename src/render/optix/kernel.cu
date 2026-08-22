#include "core/affine_transform.h"
#include "core/assert.h"
#include "core/math.h"
#include "core/matrix4x4.h"
#include "core/octahedral.h"
#include "core/vec3f.h"
#include "device/kernel_params.h"
#include "render/bsdf.h"
#include "render/camera.h"
#include "render/ray_cone.h"
#include "render/virtual_texture_interop.h"
#include <optix_device.h>
#include <optix_types.h>
#include <vector_types.h>

namespace ybi
{

// #define YBI_DEBUG_LOD_COLOR 1

struct Ray
{
    Vec3f pos;
    Vec3f dir;
};

struct Camera
{
};

struct RNG
{
    YBI_GPU float Uniform()
    {
        return 0.f;
    }
};

extern "C"
{
    __constant__ KernelParams params;
}

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
    AffineTransform worldToObject;
    uint32_t primIndex;
    uint32_t clusterId;
    uint32_t id;
    IntersectionType type;
};

YBI_GPU static Vec4f SamplePointVirtualTexture(DeviceTextureDescriptor descriptor,
                                               uint32_t pageShift,
                                               Vec2f uv,
                                               uint32_t mip,
                                               UnpackedFeedbackRequest *request = nullptr,
                                               uint32_t udimTableIndex = 0,
                                               uint32_t base = 0,
                                               uint32_t offset = 0)
{
    const uint32_t pageSize = (1u << pageShift);
    Vec2f normalizedUv = uv - Floor(uv);
    // Floating-point cancellation can make fract(uv) exactly 1 for values just below an
    // integer, particularly for negative or large UVs. Keep the repeated coordinate in
    // [0, 1) so it cannot produce a one-past-the-end texel or page.
    constexpr float maxNormalizedUv = 0x1.fffffep-1f;
    normalizedUv.x = Clamp(normalizedUv.x, 0.f, maxNormalizedUv);
    normalizedUv.y = Clamp(normalizedUv.y, 0.f, maxNormalizedUv);

    uint32_t mipWidth = descriptor.width >> mip;
    uint32_t mipHeight = descriptor.height >> mip;

    if (request)
    {
        request->textureIndex = descriptor.textureIndex;
        request->tileX = 0;
        request->tileY = 0;
        request->mip = mip;
        request->pageDescriptorIndex = ~0u;
    }

    if (mipWidth < pageSize || mipHeight < pageSize)
    {
        return descriptor.fallback;
    }

    Vec2f texel = normalizedUv * Vec2f(mipWidth, mipHeight);
    const Vec2f page = Floor(texel / float(pageSize));
    const Vec2f texelInPage = texel - page * pageSize;

    const Vec2f pageUv = texelInPage / pageSize;

    const uint32_t pageOffset = CalculatePageDescriptorOffset(descriptor.width,
                                                              descriptor.height,
                                                              mip,
                                                              pageShift,
                                                              static_cast<uint32_t>(page.x),
                                                              static_cast<uint32_t>(page.y));
    if (pageOffset >= descriptor.pageDescriptorCount)
    {
        printf("%u %u, %u: %u %u %f %f %f %f, %u %u\n",
               pageOffset,
               descriptor.pageDescriptorCount,
               udimTableIndex,
               descriptor.width,
               descriptor.height,
               page.x,
               page.y,
               uv.x,
               uv.y,
               base,
               offset);
        assert(pageOffset < descriptor.pageDescriptorCount);
    }
    const uint32_t pageIndex = descriptor.pageDescriptorOffset + pageOffset;

    if (request)
    {
        request->tileX = static_cast<uint32_t>(page.x);
        request->tileY = static_cast<uint32_t>(page.y);
        request->pageDescriptorIndex = pageIndex;
    }

    const PageDescriptor pageDescriptor = params.pageDescriptors[pageIndex];

    if (pageDescriptor.bits == ~0u)
    {
        return descriptor.fallback;
    }

    const uint32_t layer = pageDescriptor.bits >> 27u;
    const uint32_t textureIndex = pageDescriptor.bits & 0x07ffffffu;

    // DeviceTexture texture = params.texturePools[textureIndex];
    uint64_t textureObject = pageDescriptor.textureObject;
    if (textureObject == 0ull)
    {
        return descriptor.fallback;
    }

    float4 value = tex2DLayered<float4>(
        reinterpret_cast<cudaTextureObject_t>(textureObject), pageUv.x, pageUv.y, layer);

    return Vec4f(value.x, value.y, value.z, value.w);
}

YBI_GPU void WriteUniqueFeedback(UnpackedFeedbackRequest request)
{
    const uint32_t pageDescriptorIndex = request.pageDescriptorIndex;
    if (pageDescriptorIndex == ~0u)
    {
        return;
    }

    uint32_t *bits = params.pageDescriptorBits;

    const uint32_t wordIndex = pageDescriptorIndex >> 5u;
    const uint32_t bitIndex = pageDescriptorIndex & 31u;

    uint32_t old = atomicOr(&bits[wordIndex], 1u << bitIndex);
    if (((old >> bitIndex) & 1u) == 0)
    {
        uint32_t requestIndex = atomicAdd(&params.feedbackRequests[0].textureIndex, 1u);

        params.feedbackRequests[requestIndex + 1] =
            PackFeedbackRequest(request.textureIndex, request.tileX, request.tileY, request.mip);
    }
}

YBI_GPU static uint32_t CalculateMaxVirtualTextureMip(DeviceTextureDescriptor descriptor,
                                                      uint32_t pageShift)
{
    const uint32_t pageSize = 1u << pageShift;
    uint32_t maxMip = 0;

    while (maxMip < 15u)
    {
        const uint32_t nextMip = maxMip + 1u;
        if ((descriptor.width >> nextMip) < pageSize || (descriptor.height >> nextMip) < pageSize)
        {
            break;
        }

        maxMip = nextMip;
    }

    return maxMip;
}

YBI_GPU static Vec3f GetLodDebugColor(float lod, uint32_t maxMip)
{
    const float normalizedLod = maxMip ? Clamp(lod / float(maxMip), 0.f, 1.f) : 0.f;

    if (normalizedLod < 0.25f)
    {
        return Lerp(Vec3f(0.f, 0.f, 1.f), Vec3f(0.f, 1.f, 1.f), normalizedLod * 4.f);
    }
    else if (normalizedLod < 0.5f)
    {
        return Lerp(Vec3f(0.f, 1.f, 1.f), Vec3f(0.f, 1.f, 0.f), (normalizedLod - 0.25f) * 4.f);
    }
    else if (normalizedLod < 0.75f)
    {
        return Lerp(Vec3f(0.f, 1.f, 0.f), Vec3f(1.f, 1.f, 0.f), (normalizedLod - 0.5f) * 4.f);
    }

    return Lerp(Vec3f(1.f, 1.f, 0.f), Vec3f(1.f, 0.f, 0.f), (normalizedLod - 0.75f) * 4.f);
}

template <typename T>
YBI_GPU static T InterpolateTriangleValues(T values[3], Vec2f bary)
{
    return values[0] * (1.f - bary.x - bary.y) + values[1] * bary.x + values[2] * bary.y;
}

struct ShadingInfo
{
    Vec3f positions[3];
    Vec3f normals[3];
    Vec2f uvs[3];
    DeviceMaterial material;
    bool hasNormals;
    bool hasUvs;
    bool hasMaterial;
};

YBI_GPU static Vec3u GetTriangleVertexIndices(const Object &object, const IntersectionInfo &info)
{
    const uint32_t offset = object.triangleIndexOffset + 3u * info.primIndex;
    return Vec3u(static_cast<uint32_t>(params.scene.triangleIndices[offset]),
                 static_cast<uint32_t>(params.scene.triangleIndices[offset + 1u]),
                 static_cast<uint32_t>(params.scene.triangleIndices[offset + 2u]));
}

YBI_GPU static void
GetTrianglePositions(const Object &object, const Vec3u &vertexIndices, Vec3f positions[3])
{
    for (uint32_t i = 0; i < 3u; i++)
    {
        positions[i] = params.scene.positions[object.positionOffset + vertexIndices[i]];
    }
}

YBI_GPU static bool IsClusteredTriangle(const Object &object, const IntersectionInfo &info)
{
#if OPTIX_VERSION >= 90000
    return object.clusterOffset != ~0u && info.clusterId != OPTIX_CLUSTER_ID_INVALID;
#else
    return false;
#endif
}

YBI_GPU static void GetClusteredTrianglePositions(const DeviceCluster &cluster,
                                                  const uint8_t *clusterData,
                                                  const Vec3u &vertexIndices,
                                                  Vec3f positions[3])
{
    for (uint32_t i = 0; i < 3u; i++)
    {
        positions[i] = cluster.DecodePosition(clusterData, vertexIndices[i]);
    }
}

YBI_GPU static bool GetTriangleNormals(const Object &object,
                                       const IntersectionInfo &info,
                                       const Vec3u &vertexIndices,
                                       Vec3f normals[3])
{
    switch (object.normalType)
    {
        case AttributeType::FaceVarying:
        {
            for (int i = 0; i < 3; i++)
            {
                Vec3f objectNormal = DecodeOctahedral(
                    params.scene.normals[object.normalOffset + 3 * info.primIndex + i]);

                normals[i] = objectNormal;
            }
            return true;
        }
        case AttributeType::Vertex:
        {
            for (uint32_t i = 0; i < 3u; i++)
            {
                Vec3f objectNormal =
                    DecodeOctahedral(params.scene.normals[object.normalOffset + vertexIndices[i]]);
                normals[i] = objectNormal;
            }
            return true;
        }
        case AttributeType::None:
        {
            return false;
        }
        default:
            assert(0);
            return false;
    }
}

YBI_GPU static bool GetClusteredTriangleNormals(const Object &object,
                                                const IntersectionInfo &info,
                                                const DeviceCluster &cluster,
                                                const uint8_t *clusterData,
                                                const Vec3u &vertexIndices,
                                                Vec3f normals[3])
{
    switch (object.normalType)
    {
        case AttributeType::FaceVarying:
        {
            for (int i = 0; i < 3; i++)
            {
                const uint32_t normalIndex = 3u * info.primIndex + static_cast<uint32_t>(i);
                normals[i] = cluster.DecodeNormal(clusterData, normalIndex);
            }
            return true;
        }
        case AttributeType::Vertex:
        {
            for (uint32_t i = 0; i < 3u; i++)
            {
                normals[i] = cluster.DecodeNormal(clusterData, vertexIndices[i]);
            }
            return true;
        }
        case AttributeType::None:
        {
            return false;
        }
        default:
            assert(0);
            return false;
    }
}

YBI_GPU static bool GetTriangleUVs(const Object &object,
                                   const IntersectionInfo &info,
                                   const Vec3u &vertexIndices,
                                   Vec2f uvs[3])
{
    switch (object.uvType)
    {
        case AttributeType::FaceVarying:
        {
            for (int i = 0; i < 3; i++)
            {
                const uint32_t cornerIndex = 3 * info.primIndex + i;
                int uvIndex = cornerIndex;
                if (object.uvIndexOffset != ~0u)
                {
                    uvIndex = params.scene.uvIndices[object.uvIndexOffset + cornerIndex];
                }

                uvs[i] = params.scene.uvs[object.uvOffset + uvIndex];
            }
            break;
        }
        case AttributeType::Vertex:
        {
            for (uint32_t i = 0; i < 3u; i++)
            {
                int uvIndex = static_cast<int>(vertexIndices[i]);
                if (object.uvIndexOffset != ~0u)
                {
                    uvIndex = params.scene.uvIndices[object.uvIndexOffset + vertexIndices[i]];
                }

                uvs[i] = params.scene.uvs[object.uvOffset + uvIndex];
            }
            break;
        }
        case AttributeType::None:
        {
            return false;
        }
        default:
            assert(0);
    }

    return true;
}

YBI_GPU static bool GetClusteredTriangleUVs(const Object &object,
                                            const IntersectionInfo &info,
                                            const DeviceCluster &cluster,
                                            const uint8_t *clusterData,
                                            const Vec3u &vertexIndices,
                                            Vec2f uvs[3])
{
    switch (object.uvType)
    {
        case AttributeType::FaceVarying:
        {
            for (uint32_t i = 0; i < 3u; i++)
            {
                const uint32_t uvIndex = 3u * info.primIndex + i;
                uvs[i] = cluster.DecodeUV(clusterData, uvIndex);
            }
            break;
        }
        case AttributeType::Vertex:
        {
            for (uint32_t i = 0; i < 3u; i++)
            {
                uvs[i] = cluster.DecodeUV(clusterData, vertexIndices[i]);
            }
            break;
        }
        case AttributeType::None:
        {
            return false;
        }
        default:
            assert(0);
    }

    return true;
}

YBI_GPU static bool
GetTriangleMaterial(const Object &object, uint32_t triangleIndex, DeviceMaterial &material)
{
    if (object.triangleMaterialOffset != ~0u)
    {
        int materialIndex =
            params.scene.triangleMaterialIndices[object.triangleMaterialOffset + triangleIndex];
        material = params.scene.materials[materialIndex];
        return true;
    }
    else if (object.materialIndex != -1)
    {
        material = params.scene.materials[object.materialIndex];
        return true;
    }
    return false;
}

YBI_GPU static bool GetClusteredTriangleMaterial(const IntersectionInfo &info,
                                                 const DeviceCluster &cluster,
                                                 const uint8_t *clusterData,
                                                 DeviceMaterial &material)
{
    const uint32_t materialIndex = cluster.DecodeMaterialID(clusterData, info.primIndex);
    if (materialIndex == DEVICE_CLUSTER_INVALID_MATERIAL_ID)
    {
        return false;
    }

    material = params.scene.materials[materialIndex];
    return true;
}

YBI_GPU static void GetTriangleShadingInfo(const Object &object,
                                           const IntersectionInfo &info,
                                           ShadingInfo &shadingInfo)
{
    const Vec3u vertexIndices = GetTriangleVertexIndices(object, info);
    GetTrianglePositions(object, vertexIndices, shadingInfo.positions);
    shadingInfo.hasNormals = GetTriangleNormals(object, info, vertexIndices, shadingInfo.normals);
    shadingInfo.hasUvs = GetTriangleUVs(object, info, vertexIndices, shadingInfo.uvs);
    shadingInfo.hasMaterial = GetTriangleMaterial(object, info.primIndex, shadingInfo.material);
}

YBI_GPU static void GetClusteredShadingInfo(const Object &object,
                                            const IntersectionInfo &info,
                                            const DeviceCluster &cluster,
                                            const uint8_t *clusterData,
                                            const Vec3u &vertexIndices,
                                            ShadingInfo &shadingInfo)
{
    GetClusteredTrianglePositions(cluster, clusterData, vertexIndices, shadingInfo.positions);
    shadingInfo.hasNormals = GetClusteredTriangleNormals(
        object, info, cluster, clusterData, vertexIndices, shadingInfo.normals);
    shadingInfo.hasUvs = GetClusteredTriangleUVs(
        object, info, cluster, clusterData, vertexIndices, shadingInfo.uvs);
    shadingInfo.hasMaterial =
        GetClusteredTriangleMaterial(info, cluster, clusterData, shadingInfo.material);
}

YBI_GPU static void
GetShadingInfo(const Object &object, const IntersectionInfo &info, ShadingInfo &shadingInfo)
{
    shadingInfo = {};
    if (IsClusteredTriangle(object, info))
    {
        const uint32_t clusterIndex = object.clusterOffset + info.clusterId;
        const DeviceCluster &cluster = params.scene.compressedClusters[clusterIndex];
        const uint8_t *clusterData = params.scene.compressedClusterBytes + cluster.byteOffset;
        const Vec3u vertexIndices = cluster.DecodeTriangle(clusterData, info.primIndex);
        GetClusteredShadingInfo(object, info, cluster, clusterData, vertexIndices, shadingInfo);
    }
    else
    {
        GetTriangleShadingInfo(object, info, shadingInfo);
    }
}

extern "C" __global__ void __miss__kernel()
{
    uint32_t low = optixGetPayload_0();
    uint32_t high = optixGetPayload_1();
    uint64_t ptr = low | (uint64_t(high) << 32u);

    IntersectionInfo *isect = reinterpret_cast<IntersectionInfo *>(ptr);
    isect->type = INTERSECTION_TYPE_NONE;
}

extern "C" __global__ void __closesthit__kernel()
{
    uint32_t low = optixGetPayload_0();
    uint32_t high = optixGetPayload_1();
    uint64_t ptr = low | (uint64_t(high) << 32u);

    IntersectionInfo *isect = reinterpret_cast<IntersectionInfo *>(ptr);

    if (optixIsTriangleHit())
    {
        const float2 bary = optixGetTriangleBarycentrics();
        float objectToWorld[12];
        float worldToObject[12];
        optixGetObjectToWorldTransformMatrix(objectToWorld);
        optixGetWorldToObjectTransformMatrix(worldToObject);

        isect->t = optixGetRayTmax();
        isect->id = optixGetInstanceId();
        isect->primIndex = optixGetPrimitiveIndex();
        isect->clusterId = ~0u;
#if OPTIX_VERSION >= 90000
        isect->clusterId = optixGetClusterId();
#endif
        isect->bary = Vec2f(bary.x, bary.y);
        isect->objectToWorld = AffineTransform(objectToWorld);
        isect->worldToObject = AffineTransform(worldToObject);
        isect->type = INTERSECTION_TYPE_TRIANGLE;
    }
    else
    {
    }
}

extern "C" __global__ void __raygen__kernel()
{
    const uint3 launchIndex = optixGetLaunchIndex();
    const uint3 launchDimensions = optixGetLaunchDimensions();

    RNG rng;
    Vec3f pixel(float(launchIndex.x) + 0.5f, float(launchIndex.y) + 0.5f, 0.f);

    float filterU = rng.Uniform();
    float filterV = rng.Uniform();

    pixel.x += filterU - 0.5f;
    pixel.y += filterV - 0.5f;

    const int maxDepth = params.maxDepth;

    Matrix4x4 cameraFromRaster = params.cameraFromRaster;
    AffineTransform renderFromCamera = params.renderFromCamera;

    // Generate ray
    Vec4f cameraPos = cameraFromRaster * Vec4f(pixel, 1.f);
    cameraPos.xyz /= cameraPos.w;

    Vec3f rayOrigin = renderFromCamera.TransformPoint(Vec3f(0.f));
    Vec3f rayDirection = Normalize(renderFromCamera.TransformVector(cameraPos.xyz));
    float tMin = 0.f;
    float tMax = 1000000.f;
    float time = 0.f;
    uint32_t mask = 0xff;

    RayCone rayCone;
    rayCone.width = 0.f;
    rayCone.spreadAngle = atanf(2.f * tanf(0.5f * params.fovRadians) / float(params.imageHeight));

    OptixTraversableHandle bvhHandle = static_cast<OptixTraversableHandle>(params.bvhHandle);

    IntersectionInfo info;
    uint64_t ptr = reinterpret_cast<uint64_t>(&info);
    uint32_t low = ptr & 0xffffffffu;
    uint32_t high = ptr >> 32u;

    optixTrace(bvhHandle,
               make_float3(rayOrigin.x, rayOrigin.y, rayOrigin.z),
               make_float3(rayDirection.x, rayDirection.y, rayDirection.z),
               tMin,
               tMax,
               time,
               mask,
               OPTIX_RAY_FLAG_NONE,
               0,
               0,
               0,
               low,
               high);

    uchar4 color = make_uchar4(0, 0, 0, 255);
    if (info.type == INTERSECTION_TYPE_TRIANGLE)
    {
#if 1
        const Object &object = params.scene.objects[info.id];

        ShadingInfo shadingInfo;
        GetShadingInfo(object, info, shadingInfo);
        Vec3f *positions = shadingInfo.positions;
        Vec3f *normals = shadingInfo.normals;

        if (!shadingInfo.hasNormals)
        {
            Vec3f n = Normalize(Cross(positions[1] - positions[0], positions[2] - positions[0]));

            normals[0] = n;
            normals[1] = n;
            normals[2] = n;
        }

        Vec2f *uvs = shadingInfo.uvs;
        Vec2f uv = InterpolateTriangleValues(uvs, info.bary);

        DeviceMaterial &material = shadingInfo.material;

        if (shadingInfo.hasMaterial && material.valid)
        {
            if (shadingInfo.hasUvs)
            {
                // Virtual texture path
                Vec2f flooredUv = Floor(uv);
                uint32_t udimTableIndex = 0u;
                if (material.udimBase != ~0u)
                {
                    const uint32_t udim = 1001 + static_cast<uint32_t>(flooredUv.x) +
                                          10 * static_cast<uint32_t>(flooredUv.y);
                    udimTableIndex = material.udimTableOffset + udim - material.udimBase;
                }
                else
                {
                    udimTableIndex = material.udimTableOffset;
                }

                DeviceTextureDescriptor textureDescriptor =
                    params.textureDescriptors[udimTableIndex];

                const uint32_t pageShift = 7;
                const Vec3f hitPosition = rayOrigin + rayDirection * info.t;
                const Vec3f cameraPosition = params.cameraPosition;

                // Calculate LOD
                for (int i = 0; i < 3; i++)
                {
                    positions[i] = info.objectToWorld.TransformPoint(positions[i]);
                }

                Vec3f r0 = info.objectToWorld.GetRow(0).xyz;
                Vec3f r1 = info.objectToWorld.GetRow(1).xyz;
                Vec3f r2 = info.objectToWorld.GetRow(2).xyz;

                Vec3f adjugate0 = Cross(r1, r2);
                Vec3f adjugate1 = Cross(r2, r0);
                Vec3f adjugate2 = Cross(r0, r1);

                Vec3f shadingNormal = Normalize(InterpolateTriangleValues(normals, info.bary));
                shadingNormal = Normalize(Vec3f(Dot(adjugate0, shadingNormal),
                                                Dot(adjugate1, shadingNormal),
                                                Dot(adjugate2, shadingNormal)));

                for (int i = 0; i < 3; i++)
                {
                    normals[i] = Normalize(Vec3f(Dot(adjugate0, normals[i]),
                                                 Dot(adjugate1, normals[i]),
                                                 Dot(adjugate2, normals[i])));
                }

                float surfaceSpreadAngle =
                    CalculatePrimaryHitUnifiedSurfaceSpreadAngle(rayCone,
                                                                 rayDirection,
                                                                 shadingNormal,
                                                                 positions[0],
                                                                 positions[1],
                                                                 positions[2],
                                                                 normals[0],
                                                                 normals[1],
                                                                 normals[2]);
                PropagateRayCone(rayCone, surfaceSpreadAngle, info.t);
                float lod = ComputeTextureLOD(positions[0],
                                              positions[1],
                                              positions[2],
                                              uvs[0],
                                              uvs[1],
                                              uvs[2],
                                              rayDirection,
                                              shadingNormal,
                                              rayCone.width,
                                              textureDescriptor.width,
                                              textureDescriptor.height);

                const uint32_t maxMip =
                    CalculateMaxVirtualTextureMip(textureDescriptor, pageShift);
                lod = Clamp(lod, 0.f, float(maxMip));

#if YBI_DEBUG_LOD_COLOR
                Vec3f value = GetLodDebugColor(lod, maxMip);
                color = make_uchar4(static_cast<unsigned char>(Clamp(value.x, 0.f, 1.f) * 255.f),
                                    static_cast<unsigned char>(Clamp(value.y, 0.f, 1.f) * 255.f),
                                    static_cast<unsigned char>(Clamp(value.z, 0.f, 1.f) * 255.f),
                                    255);
#else
                uint32_t baseMip = uint32_t(Floor(lod));
                uint32_t nextMip = Min(baseMip + 1, maxMip);

                // TODO: don't hardcode
                UnpackedFeedbackRequest lowR, highR;
                Vec4f lowerValue = SamplePointVirtualTexture(textureDescriptor,
                                                             pageShift,
                                                             uv,
                                                             baseMip,
                                                             &lowR,
                                                             udimTableIndex,
                                                             material.udimBase,
                                                             material.udimTableOffset);
                Vec4f higherValue = SamplePointVirtualTexture(textureDescriptor,
                                                              pageShift,
                                                              uv,
                                                              nextMip,
                                                              &highR,
                                                              udimTableIndex,
                                                              material.udimBase,
                                                              material.udimTableOffset);

                WriteUniqueFeedback(lowR);
                WriteUniqueFeedback(highR);

                float t = lod - float(baseMip);
                Vec4f value = Lerp(lowerValue, higherValue, t);
                color = make_uchar4(static_cast<unsigned char>(Clamp(value.x, 0.f, 1.f) * 255.f),
                                    static_cast<unsigned char>(Clamp(value.y, 0.f, 1.f) * 255.f),
                                    static_cast<unsigned char>(Clamp(value.z, 0.f, 1.f) * 255.f),
                                    static_cast<unsigned char>(Clamp(value.w, 0.f, 1.f) * 255.f));
#endif
            }
        }
#else
        uint32_t hash = info.clusterId;
        hash ^= hash >> 16u;
        hash *= 0x7feb352du;
        hash ^= hash >> 15u;
        hash *= 0x846ca68bu;
        hash ^= hash >> 16u;

        color.x = static_cast<unsigned char>(64u + (hash & 0xbfu));
        color.y = static_cast<unsigned char>(64u + ((hash >> 8u) & 0xbfu));
        color.z = static_cast<unsigned char>(64u + ((hash >> 16u) & 0xbfu));
#endif
    }

    uchar4 *image = reinterpret_cast<uchar4 *>(params.imageHandle);
    const uint32_t pixelIndex = launchIndex.x + launchIndex.y * launchDimensions.x;
    image[pixelIndex] = color;
}

} // namespace ybi
