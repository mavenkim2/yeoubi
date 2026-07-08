#include "core/affine_transform.h"
#include "core/assert.h"
#include "core/math.h"
#include "core/matrix4x4.h"
#include "core/vec3f.h"
#include "device/kernel_params.h"
#include "render/ray_cone.h"
#include "render/virtual_texture_interop.h"
#include <optix_device.h>
#include <optix_types.h>
#include <vector_types.h>

namespace ybi
{

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
    uint32_t primIndex;
    uint32_t id;
    IntersectionType type;
};

YBI_GPU static Vec4f SamplePointVirtualTexture(DeviceTextureDescriptor descriptor,
                                               uint32_t pageShift,
                                               Vec2f uv,
                                               uint32_t mip,
                                               UnpackedFeedbackRequest *request = nullptr)
{
    const uint32_t pageSize = (1u << pageShift);
    const Vec2f normalizedUv = uv - Floor(uv);

    const uint32_t mipWidth = descriptor.width >> mip;
    const uint32_t mipHeight = descriptor.height >> mip;

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
    const Vec2f texelInPage = (texel + 0.5f) - page * pageSize;

    const Vec2f pageUv = texelInPage / pageSize;

    const uint32_t pageOffset = CalculatePageDescriptorOffset(descriptor.width,
                                                              descriptor.height,
                                                              mip,
                                                              pageShift,
                                                              static_cast<uint32_t>(page.x),
                                                              static_cast<uint32_t>(page.y));
    assert(pageOffset < descriptor.pageDescriptorCount);
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

template <typename T>
YBI_GPU static T InterpolateTriangleValues(T values[3], Vec2f bary)
{
    return values[0] * (1.f - bary.x - bary.y) + values[1] * bary.x + values[2] * bary.y;
}

YBI_GPU static void
GetTrianglePositions(const Object &object, const IntersectionInfo &info, Vec3f positions[3])
{
    const uint32_t triangleIndexOffset = object.triangleIndexOffset + 3 * info.primIndex;
    for (int i = 0; i < 3; i++)
    {
        const int vertexIndex = params.scene.triangleIndices[triangleIndexOffset + i];
        positions[i] = params.scene.positions[object.positionOffset + vertexIndex];
    }
}

YBI_GPU static bool
GetTriangleNormals(const Object &object, const IntersectionInfo &info, Vec3f normals[3])
{
    if (object.normalOffset == ~0u || params.scene.normals == nullptr)
    {
        return false;
    }

    switch (object.normalType)
    {
        case AttributeType::FaceVarying:
        {
            for (int i = 0; i < 3; i++)
            {
                normals[i] = params.scene.normals[object.normalOffset + 3 * info.primIndex + i];
            }
            return true;
        }
        case AttributeType::Vertex:
        {
            const uint32_t triangleIndexOffset = object.triangleIndexOffset + 3 * info.primIndex;
            for (int i = 0; i < 3; i++)
            {
                const int vertexIndex = params.scene.triangleIndices[triangleIndexOffset + i];
                normals[i] = params.scene.normals[object.normalOffset + vertexIndex];
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

YBI_GPU static bool
GetTriangleUVs(const Object &object, const IntersectionInfo &info, Vec2f uvs[3])
{
    switch (object.uvType)
    {
        case AttributeType::FaceVarying:
        {
            for (int i = 0; i < 3; i++)
            {
                uvs[i] = params.scene.float2s[object.uvOffset + 3 * info.primIndex + i];
            }
            break;
        }
        case AttributeType::Vertex:
        {
            const uint32_t triangleIndexOffset = object.triangleIndexOffset + 3 * info.primIndex;
            for (int i = 0; i < 3; i++)
            {
                const int vertexIndex = params.scene.triangleIndices[triangleIndexOffset + i];
                uvs[i] = params.scene.float2s[object.uvOffset + vertexIndex];
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
GetTriangleMaterial(const Object &object, const IntersectionInfo &info, DeviceMaterial &material)
{
    if (object.triangleMaterialOffset != ~0u)
    {
        int materialIndex =
            params.scene.triangleMaterialIndices[object.triangleMaterialOffset + info.primIndex];
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
        isect->t = optixGetRayTmax();
        isect->id = optixGetInstanceId();
        isect->primIndex = optixGetPrimitiveIndex();
        isect->bary = Vec2f(bary.x, bary.y);
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
        const Object &object = params.scene.objects[info.id];

        Vec3f positions[3];
        Vec3f normals[3];

        GetTrianglePositions(object, info, positions);

        bool hasNormals = GetTriangleNormals(object, info, normals);

        if (!hasNormals)
        {
            Vec3f n = Normalize(Cross(positions[1] - positions[0], positions[2] - positions[0]));

            normals[0] = n;
            normals[1] = n;
            normals[2] = n;
        }

        Vec2f uvs[3];
        bool hasUvs = GetTriangleUVs(object, info, uvs);
        Vec2f uv = InterpolateTriangleValues(uvs, info.bary);

        DeviceMaterial material;
        bool hasMaterial = GetTriangleMaterial(object, info, material);

        if (hasMaterial)
        {
            if (hasUvs)
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
                Vec3f shadingNormal = InterpolateTriangleValues(normals, info.bary);

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
                                              params.imageWidth,
                                              params.imageHeight);

                const uint32_t maxMip =
                    CalculateMaxVirtualTextureMip(textureDescriptor, pageShift);
                lod = Clamp(lod, 0.f, float(maxMip));

                uint32_t baseMip = uint32_t(Floor(lod));
                uint32_t nextMip = Min(baseMip + 1, maxMip);

#if 1
                // TODO: don't hardcode
                UnpackedFeedbackRequest lowR, highR;
                Vec4f lowerValue =
                    SamplePointVirtualTexture(textureDescriptor, pageShift, uv, baseMip, &lowR);
                Vec4f higherValue =
                    SamplePointVirtualTexture(textureDescriptor, pageShift, uv, nextMip, &highR);

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
    }

    uchar4 *image = reinterpret_cast<uchar4 *>(params.imageHandle);
    const uint32_t pixelIndex = launchIndex.x + launchIndex.y * launchDimensions.x;
#if 0
    if (info.type == INTERSECTION_TYPE_TRIANGLE)
    {
        uint32_t hash = info.id;
        hash ^= hash >> 16;
        hash *= 0x7feb352du;
        hash ^= hash >> 15;
        hash *= 0x846ca68bu;
        hash ^= hash >> 16;

        color.x = static_cast<unsigned char>(64u + (hash & 0xbfu));
        color.y = static_cast<unsigned char>(64u + ((hash >> 8u) & 0xbfu));
        color.z = static_cast<unsigned char>(64u + ((hash >> 16u) & 0xbfu));
    }
#endif
    image[pixelIndex] = color;

    // for (int depth = 0; depth < maxDepth; depth++)
    // {
    // optixTraverse(, );
    // }

    // Generate a ray from the camera
    // Intersect the scene using the ray
}

} // namespace ybi
