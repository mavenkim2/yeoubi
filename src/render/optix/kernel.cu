#include "core/affine_transform.h"
#include "core/assert.h"
#include "core/math.h"
#include "core/matrix4x4.h"
#include "core/octahedral.h"
#include "core/vec3f.h"
#include "device/kernel_params.h"
#include "render/bsdf/openpbr/glossy_diffuse.h"
#include "render/camera.h"
#include "render/device_attributes.h"
#include "render/ray_cone.h"
#include "render/ray_cone_utils.h"
#include "render/rng.h"
#include "render/sampling.h"
#include "render/shaders/material_shaders.h"
#include "render/shaders/shading.h"
#include "render/virtual_texture_interop.h"
#include "render/virtual_texture_sampling.h"
#include "scene/openpbr_parameters.h"
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
GetTriangleMaterial(const Object &object, uint32_t triangleIndex, uint32_t &materialIndex)
{
    if (object.triangleMaterialOffset != ~0u)
    {
        materialIndex =
            params.scene.triangleMaterialIndices[object.triangleMaterialOffset + triangleIndex];
        return true;
    }
    else if (object.materialIndex != -1)
    {
        materialIndex = object.materialIndex;
        return true;
    }
    return false;
}

YBI_GPU static bool GetClusteredTriangleMaterial(const IntersectionInfo &info,
                                                 const DeviceCluster &cluster,
                                                 const uint8_t *clusterData,
                                                 uint32_t &materialIndex)
{
    materialIndex = cluster.DecodeMaterialID(clusterData, info.primIndex);
    if (materialIndex == DEVICE_CLUSTER_INVALID_MATERIAL_ID)
    {
        return false;
    }

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
    shadingInfo.hasMaterial =
        GetTriangleMaterial(object, info.primIndex, shadingInfo.materialIndex);
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
        GetClusteredTriangleMaterial(info, cluster, clusterData, shadingInfo.materialIndex);
}

YBI_GPU static void
GetShadingInfo(const Object &object, const IntersectionInfo &info, ShadingInfo &shadingInfo)
{
    shadingInfo = {};
    if (IsClusteredTriangle(object.clusterOffset, info.clusterId))
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

    shadingInfo.bary = info.bary;
    CalculateTriangleHitInfo(shadingInfo);
    shadingInfo.objectId = info.id;
    shadingInfo.clusterId = info.clusterId;
    shadingInfo.primIndex = info.primIndex;
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

    RNG rng = RNG::Init(Vec2u(launchIndex.x, launchIndex.y), 0u);
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

        if (shadingInfo.hasMaterial)
        {
            const MaterialHandle &handle = params.scene.materialHandles[shadingInfo.materialIndex];

            if (handle.isGraph)
            {
                const DeviceMaterial2 &material = params.scene.deviceMaterial2s[handle.index];
                const Vec3f hitPosition = rayOrigin + rayDirection * info.t;
                const Vec3f cameraPosition = params.cameraPosition;

                shadingInfo.TransformToWorld(info.objectToWorld);

                // Calculate LOD

                float surfaceSpreadAngle =
                    CalculatePrimaryHitUnifiedSurfaceSpreadAngle(rayCone,
                                                                 rayDirection,
                                                                 shadingInfo.shadingNormal,
                                                                 positions[0],
                                                                 positions[1],
                                                                 positions[2],
                                                                 normals[0],
                                                                 normals[1],
                                                                 normals[2]);
                PropagateRayCone(rayCone, surfaceSpreadAngle, info.t);
                shadingInfo.rayDirection = rayDirection;
                shadingInfo.rayCone = rayCone;

                const Vec3f baseColor =
                    EvaluateShader<Vec3f>(params, shadingInfo, material.base_color);

                // Vec3f baseColor(0.f);
                color =
                    make_uchar4(static_cast<unsigned char>(Clamp(baseColor.x, 0.f, 1.f) * 255.f),
                                static_cast<unsigned char>(Clamp(baseColor.y, 0.f, 1.f) * 255.f),
                                static_cast<unsigned char>(Clamp(baseColor.z, 0.f, 1.f) * 255.f),
                                255);
            }
            // TODO: get rid of this legacy path
            else
            {
                DeviceMaterial &material = params.scene.materials[handle.index];
                if (material.valid && shadingInfo.hasUvs)
                {
                    shadingInfo.TransformToWorld(info.objectToWorld);

                    float surfaceSpreadAngle =
                        CalculatePrimaryHitUnifiedSurfaceSpreadAngle(rayCone,
                                                                     rayDirection,
                                                                     shadingInfo.shadingNormal,
                                                                     positions[0],
                                                                     positions[1],
                                                                     positions[2],
                                                                     normals[0],
                                                                     normals[1],
                                                                     normals[2]);
                    PropagateRayCone(rayCone, surfaceSpreadAngle, info.t);

                    shadingInfo.rayDirection = rayDirection;
                    shadingInfo.rayCone = rayCone;

                    const MaterialTextureHandle textureHandle = {material.udimBase,
                                                                 material.udimTableOffset};
                    const Vec4f value =
                        VirtualTextureSampleImage(params, textureHandle, shadingInfo);
                    color =
                        make_uchar4(static_cast<unsigned char>(Clamp(value.x, 0.f, 1.f) * 255.f),
                                    static_cast<unsigned char>(Clamp(value.y, 0.f, 1.f) * 255.f),
                                    static_cast<unsigned char>(Clamp(value.z, 0.f, 1.f) * 255.f),
                                    static_cast<unsigned char>(Clamp(value.w, 0.f, 1.f) * 255.f));
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
        }
    }

    uchar4 *image = reinterpret_cast<uchar4 *>(params.imageHandle);
    const uint32_t pixelIndex = launchIndex.x + launchIndex.y * launchDimensions.x;
    image[pixelIndex] = color;
}

#if 1
extern "C" __global__ void __raygen__unidirectional()
{
    const uint3 launchIndex = optixGetLaunchIndex();
    const uint3 launchDimensions = optixGetLaunchDimensions();

    RNG rng = RNG::Init(Vec2u(launchIndex.x, launchIndex.y), 0u);
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

    Vec3f environment(0.5f);
    Vec3f radiance(0.f);
    Vec3f throughput(1.f);

    for (uint32_t depth = 0; depth < params.maxDepth + 1; depth++)
    {
        float time = 0.f;

        IntersectionInfo info;
        uint64_t ptr = reinterpret_cast<uint64_t>(&info);
        uint32_t low = ptr & 0xffffffffu;
        uint32_t high = ptr >> 32u;

        // Interect ray against scene
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

        if (info.type == INTERSECTION_TYPE_NONE)
        {
            radiance += throughput * environment;
            break;
        }

        if (depth == params.maxDepth)
        {
            break;
        }

        // Sample BSDF
        // NOTE: includes PDF, cosine term
        const Object &object = params.scene.objects[info.id];

        ShadingInfo shadingInfo;
        GetShadingInfo(object, info, shadingInfo);

        Vec2f bsdfU = rng.Uniform2D();

        if (shadingInfo.hasMaterial)
        {
            const MaterialHandle &handle = params.scene.materialHandles[shadingInfo.materialIndex];
            YBI_ASSERT(handle.isGraph);

            const DeviceMaterial2 &material = params.scene.deviceMaterial2s[handle.index];
            const Vec3f hitPosition = rayOrigin + rayDirection * info.t;
            const Vec3f cameraPosition = params.cameraPosition;

            shadingInfo.TransformToWorld(info.objectToWorld);

            UpdateRayCone(rayCone, rayDirection, info.t, shadingInfo, depth);
            shadingInfo.rayDirection = rayDirection;
            shadingInfo.rayCone = rayCone;

            // Sample BSDF
            Vec3f direction = SampleCosineHemisphere(bsdfU);
            const Vec3f baseColor =
                EvaluateShader<Vec3f>(params, shadingInfo, material.base_color);

            throughput = throughput * baseColor;

            // Transform ray back to world space
        }
        else
        {
            // Assume no-op, change ray origin and preserve ray direction
        }

        if (throughput.x == 0.f && throughput.y == 0.f && throughput.z == 0.f)
        {
            break;
        }
    }
}
#endif

} // namespace ybi
