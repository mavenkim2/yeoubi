#include "core/affine_transform.h"
#include "core/matrix4x4.h"
#include "core/vec3f.h"
#include "device/kernel_params.h"
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
    float u;
    float v;
    uint32_t primIndex;
    uint32_t id;
    IntersectionType type;
};

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
        isect->u = bary.x;
        isect->v = bary.y;
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

    uchar4 *image = reinterpret_cast<uchar4 *>(params.imageHandle);
    const uint32_t pixelIndex = launchIndex.x + launchIndex.y * launchDimensions.x;
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
    image[pixelIndex] = color;

    // for (int depth = 0; depth < maxDepth; depth++)
    // {
    // optixTraverse(, );
    // }

    // Generate a ray from the camera
    // Intersect the scene using the ray
}

} // namespace ybi
