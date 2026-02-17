#include <optix.h>
#include <optix_device.h>

struct LaunchParams
{
    struct WireframeConfig
    {
        float lineWidth;
        float lineFeather;
        float edgeDarkness;
        float padding;
    };

    OptixTraversableHandle traversable;
    unsigned long long image;
    int width;
    int height;
    float3 cameraOrigin;
    float3 cameraU;
    float3 cameraV;
    float3 cameraW;
    WireframeConfig wireframe;
};

extern "C"
{
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

static __forceinline__ __device__ float3 SkyColor(const float3 &direction)
{
    const float t = 0.5f * (direction.y + 1.0f);
    const float3 top = make_float3(0.7f, 0.8f, 1.0f);
    const float3 bottom = make_float3(0.2f, 0.25f, 0.35f);
    return make_float3((1.0f - t) * top.x + t * bottom.x,
                       (1.0f - t) * top.y + t * bottom.y,
                       (1.0f - t) * top.z + t * bottom.z);
}

// Kept for debug parity with previous shader mode; intentionally unused.
static __forceinline__ __device__ float3 ShadeBarycentric(const float3 &barycentrics,
                                                          unsigned int primitive)
{
    const float hash = (float)((primitive * 1664525u + 1013904223u) & 255u) / 255.0f;
    return make_float3(0.25f + 0.6f * barycentrics.x + 0.15f * hash,
                       0.25f + 0.6f * barycentrics.y + 0.1f * (1.0f - hash),
                       0.25f + 0.6f * barycentrics.z);
}

static __forceinline__ __device__ float SmoothStep(float edge0, float edge1, float x)
{
    const float denom = fmaxf(edge1 - edge0, 1e-6f);
    const float t = fminf(fmaxf((x - edge0) / denom, 0.0f), 1.0f);
    return t * t * (3.0f - 2.0f * t);
}

static __forceinline__ __device__ float3 ComputeDirection(const uint3 &launchIndex,
                                                          const uint3 &launchDims,
                                                          const float2 &pixelOffset)
{
    const float2 ndc =
        make_float2(((float)launchIndex.x + pixelOffset.x) / (float)launchDims.x * 2.0f - 1.0f,
                    ((float)launchIndex.y + pixelOffset.y) / (float)launchDims.y * 2.0f - 1.0f);

    return Normalize3(
        make_float3(params.cameraU.x * ndc.x + params.cameraV.x * ndc.y + params.cameraW.x,
                    params.cameraU.y * ndc.x + params.cameraV.y * ndc.y + params.cameraW.y,
                    params.cameraU.z * ndc.x + params.cameraV.z * ndc.y + params.cameraW.z));
}

static __forceinline__ __device__ unsigned int TraceColor(const float3 &origin,
                                                          const float3 &direction)
{
    unsigned int packedColor = 0;
    optixTrace(params.traversable,
               origin,
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
    return packedColor;
}

extern "C" __global__ void __raygen__primary()
{
    const uint3 launchIndex = optixGetLaunchIndex();
    const uint3 launchDims = optixGetLaunchDimensions();
    const unsigned int pixelIndex = launchIndex.y * launchDims.x + launchIndex.x;

    const float2 centerOffset = make_float2(0.5f, 0.5f);
    const float3 centerDirection = ComputeDirection(launchIndex, launchDims, centerOffset);
    const unsigned int packedColor = TraceColor(params.cameraOrigin, centerDirection);
    uchar4 *image = reinterpret_cast<uchar4 *>(params.image);
    image[pixelIndex] = make_uchar4((unsigned char)(packedColor & 255u),
                                    (unsigned char)((packedColor >> 8) & 255u),
                                    (unsigned char)((packedColor >> 16) & 255u),
                                    255u);
}

extern "C" __global__ void __miss__primary()
{
    const float3 direction = optixGetWorldRayDirection();
    const float3 color = SkyColor(direction);
    optixSetPayload_0(PackColor(color.x, color.y, color.z));
}

extern "C" __global__ void __anyhit__primary() {}

extern "C" __global__ void __closesthit__primary()
{
    float3 color = make_float3(0.86f, 0.86f, 0.86f);
    const unsigned int hitKind = optixGetHitKind();
    if (hitKind == OPTIX_HIT_KIND_TRIANGLE_FRONT_FACE ||
        hitKind == OPTIX_HIT_KIND_TRIANGLE_BACK_FACE)
    {
        const float2 bary = optixGetTriangleBarycentrics();
        const float3 barycentrics = make_float3(1.0f - bary.x - bary.y, bary.x, bary.y);
        const float edgeDistance = fminf(barycentrics.x, fminf(barycentrics.y, barycentrics.z));
        const float edgeAlpha =
            1.0f - SmoothStep(params.wireframe.lineWidth,
                              params.wireframe.lineWidth + params.wireframe.lineFeather,
                              edgeDistance);
        const float edgeValue = fminf(fmaxf(params.wireframe.edgeDarkness, 0.0f), 1.0f);
        const float3 edgeColor = make_float3(edgeValue, edgeValue, edgeValue);
        color = make_float3(color.x * (1.0f - edgeAlpha) + edgeColor.x * edgeAlpha,
                            color.y * (1.0f - edgeAlpha) + edgeColor.y * edgeAlpha,
                            color.z * (1.0f - edgeAlpha) + edgeColor.z * edgeAlpha);
        (void)ShadeBarycentric(barycentrics, optixGetPrimitiveIndex());
    }
    else
    {
        color = make_float3(0.95f, 0.90f, 0.72f);
    }
    optixSetPayload_0(PackColor(color.x, color.y, color.z));
}
