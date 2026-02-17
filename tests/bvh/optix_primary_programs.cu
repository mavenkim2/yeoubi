#include <optix.h>
#include <optix_device.h>

struct LaunchParams
{
    OptixTraversableHandle traversable;
    unsigned long long image;
    int width;
    int height;
    float3 cameraOrigin;
    float3 cameraU;
    float3 cameraV;
    float3 cameraW;
};

extern "C" {
__constant__ LaunchParams params;
}

static __forceinline__ __device__ float3 Normalize3(const float3 &v)
{
    const float invLen = rsqrtf(v.x * v.x + v.y * v.y + v.z * v.z + 1e-20f);
    return make_float3(v.x * invLen, v.y * invLen, v.z * invLen);
}

static __forceinline__ __device__ unsigned int PackColor(float r, float g, float b)
{
    const unsigned int ru = (unsigned int)(fminf(fmaxf(r, 0.0f), 1.0f) * 255.0f + 0.5f);
    const unsigned int gu = (unsigned int)(fminf(fmaxf(g, 0.0f), 1.0f) * 255.0f + 0.5f);
    const unsigned int bu = (unsigned int)(fminf(fmaxf(b, 0.0f), 1.0f) * 255.0f + 0.5f);
    return ru | (gu << 8) | (bu << 16);
}

extern "C" __global__ void __raygen__primary()
{
    const uint3 launchIndex = optixGetLaunchIndex();
    const uint3 launchDims = optixGetLaunchDimensions();
    const unsigned int pixelIndex = launchIndex.y * launchDims.x + launchIndex.x;

    const float2 ndc = make_float2(
        ((float)launchIndex.x + 0.5f) / (float)launchDims.x * 2.0f - 1.0f,
        ((float)launchIndex.y + 0.5f) / (float)launchDims.y * 2.0f - 1.0f);

    const float3 direction = Normalize3(make_float3(params.cameraU.x * ndc.x + params.cameraV.x * ndc.y +
                                                         params.cameraW.x,
                                                     params.cameraU.y * ndc.x + params.cameraV.y * ndc.y +
                                                         params.cameraW.y,
                                                     params.cameraU.z * ndc.x + params.cameraV.z * ndc.y +
                                                         params.cameraW.z));

    unsigned int packedColor = 0;
    optixTrace(params.traversable,
               params.cameraOrigin,
               direction,
               0.001f,
               1e20f,
               0.0f,
               OptixVisibilityMask(255),
               OPTIX_RAY_FLAG_NONE,
               0,
               1,
               0,
               packedColor);

    uchar4 *image = reinterpret_cast<uchar4 *>(params.image);
    image[pixelIndex] = make_uchar4((unsigned char)(packedColor & 255u),
                                    (unsigned char)((packedColor >> 8) & 255u),
                                    (unsigned char)((packedColor >> 16) & 255u),
                                    255u);
}

extern "C" __global__ void __miss__primary()
{
    const float3 direction = optixGetWorldRayDirection();
    const float t = 0.5f * (direction.y + 1.0f);
    const float3 top = make_float3(0.7f, 0.8f, 1.0f);
    const float3 bottom = make_float3(0.2f, 0.25f, 0.35f);
    const float3 color = make_float3((1.0f - t) * top.x + t * bottom.x,
                                     (1.0f - t) * top.y + t * bottom.y,
                                     (1.0f - t) * top.z + t * bottom.z);
    optixSetPayload_0(PackColor(color.x, color.y, color.z));
}

extern "C" __global__ void __closesthit__primary()
{
    const float2 bary = optixGetTriangleBarycentrics();
    const float b0 = 1.0f - bary.x - bary.y;
    const float b1 = bary.x;
    const float b2 = bary.y;
    const unsigned int primitive = optixGetPrimitiveIndex();
    const float hash = (float)((primitive * 1664525u + 1013904223u) & 255u) / 255.0f;
    const float3 color = make_float3(0.25f + 0.6f * b0 + 0.15f * hash,
                                     0.25f + 0.6f * b1 + 0.1f * (1.0f - hash),
                                     0.25f + 0.6f * b2);
    optixSetPayload_0(PackColor(color.x, color.y, color.z));
}
