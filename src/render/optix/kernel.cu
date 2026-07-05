#include "core/affine_transform.h"
#include "core/assert.h"
#include "core/math.h"
#include "core/matrix4x4.h"
#include "core/vec3f.h"
#include "device/kernel_params.h"
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
                                               uint32_t mip)
{
    const uint32_t pageSize = (1u << pageShift);
    const Vec2f normalizedUv = uv - Floor(uv);

    const uint32_t mipWidth = descriptor.width >> mip;
    const uint32_t mipHeight = descriptor.height >> mip;

    if (mipWidth < pageSize || mipHeight < pageSize)
    {
        return descriptor.fallback;
    }

    Vec2f texel = normalizedUv * Vec2f(mipWidth, mipHeight);
    const Vec2f page = Floor(texel / float(pageSize));
    const Vec2f texelInPage = (texel + 0.5f) - page * pageSize;

    const Vec2f pageUv = texelInPage / pageSize;

    const uint32_t pageIndex = descriptor.pageDescriptorOffset +
                               CalculatePageDescriptorOffset(descriptor.width,
                                                             descriptor.height,
                                                             mip,
                                                             pageShift,
                                                             static_cast<uint32_t>(page.x),
                                                             static_cast<uint32_t>(page.y));
    const PageDescriptor pageDescriptor = params.pageDescriptors[pageIndex];

    if (pageDescriptor.bits == ~0u)
    {
        return descriptor.fallback;
    }

    const uint32_t layer = pageDescriptor.bits >> 27u;
    const uint32_t textureIndex = pageDescriptor.bits & 0x07ffffffu;

    DeviceTexture texture = params.texturePools[textureIndex];

    float4 value = tex2DLayered<float4>(
        reinterpret_cast<cudaTextureObject_t>(texture.textureObject), pageUv.x, pageUv.y, layer);
    return Vec4f(value.x, value.y, value.z, value.w);
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

        Vec2f uvs[3];
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
                // TODO: need triangle indices
                // for (int i = 0; i < 3; i++)
                // {
                // }
                // break;
            }
            default:
                YBI_ASSERT(0);
        }

        Vec2f uv = uvs[0] * (1.f - info.bary.x - info.bary.y) + uvs[1] * info.bary.x +
                   uvs[2] * info.bary.y;

        // Virtual texture path
        Vec2f flooredUv = Floor(uv);
        const uint32_t udim =
            1001 + static_cast<uint32_t>(flooredUv.x) + 10 * static_cast<uint32_t>(flooredUv.y);
        const uint32_t udimTableIndex = object.udimTableOffset + udim - object.udimBase;

        DeviceTextureDescriptor textureDescriptor = params.textureDescriptors[udimTableIndex];
        float lod = 0.f;
        uint32_t baseMip = uint32_t(Floor(lod));
        uint32_t nextMip = baseMip + 1;

        // TODO: don't hardcode
        const uint32_t pageShift = 7;
        Vec4f lowerValue = SamplePointVirtualTexture(textureDescriptor, pageShift, uv, baseMip);
        Vec4f higherValue = SamplePointVirtualTexture(textureDescriptor, pageShift, uv, nextMip);

        float t = lod - float(baseMip);
        Vec4f value = Lerp(lowerValue, higherValue, t);
        color = make_uchar4(static_cast<unsigned char>(Clamp(value.x, 0.f, 1.f) * 255.f),
                            static_cast<unsigned char>(Clamp(value.y, 0.f, 1.f) * 255.f),
                            static_cast<unsigned char>(Clamp(value.z, 0.f, 1.f) * 255.f),
                            static_cast<unsigned char>(Clamp(value.w, 0.f, 1.f) * 255.f));
    }

    uchar4 *image = reinterpret_cast<uchar4 *>(params.imageHandle);
    const uint32_t pixelIndex = launchIndex.x + launchIndex.y * launchDimensions.x;
#if 0
    uchar4 color = make_uchar4(0, 0, 0, 255);
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
