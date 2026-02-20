#include <optix.h>
#include <optix_device.h>

struct LaunchParams
{
    struct InstanceGeomRef
    {
        unsigned long long positions;
        unsigned long long indices;
        int numPositions;
        int numIndices;
    };

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
    int integrator;
    int spp;
    float aoBias;
    float aoMaxDistance;
    unsigned long long instanceGeomRefs;
    int instanceGeomRefCount;
};

struct HitgroupData
{
    unsigned long long positions;
    unsigned long long indices;
    int numPositions;
    int numIndices;
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

static __forceinline__ __device__ float3 Cross3(const float3 &a, const float3 &b)
{
    return make_float3(
        a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
}

static __forceinline__ __device__ unsigned int Hash32(unsigned int x)
{
    x ^= x >> 16;
    x *= 0x7feb352du;
    x ^= x >> 15;
    x *= 0x846ca68bu;
    x ^= x >> 16;
    return x;
}

static __forceinline__ __device__ float Random01(unsigned int &state)
{
    state = Hash32(state + 0x9e3779b9u);
    return float(state & 0x00ffffffu) / float(0x01000000u);
}

static __forceinline__ __device__ void BuildOrthonormalBasis(const float3 &n, float3 &t, float3 &b)
{
    const float3 up = fabsf(n.z) < 0.999f ? make_float3(0.0f, 0.0f, 1.0f) : make_float3(0.0f, 1.0f, 0.0f);
    t = Normalize3(Cross3(up, n));
    b = Normalize3(Cross3(n, t));
}

static __forceinline__ __device__ float3 SampleCosineHemisphere(float u1, float u2)
{
    const float r = sqrtf(fmaxf(u1, 0.0f));
    const float phi = 6.28318530718f * u2;
    const float x = r * cosf(phi);
    const float y = r * sinf(phi);
    const float z = sqrtf(fmaxf(1.0f - u1, 0.0f));
    return make_float3(x, y, z);
}

static __forceinline__ __device__ float3 FaceForward(const float3 &normal, const float3 &referenceDirection)
{
    const float d = normal.x * referenceDirection.x + normal.y * referenceDirection.y + normal.z * referenceDirection.z;
    if (d < 0.0f)
    {
        return normal;
    }
    return make_float3(-normal.x, -normal.y, -normal.z);
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
                    1.0f - ((float)launchIndex.y + pixelOffset.y) / (float)launchDims.y * 2.0f);

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

static __forceinline__ __device__ unsigned int TraceOcclusion(const float3 &origin,
                                                              const float3 &direction,
                                                              float tMin,
                                                              float tMax)
{
    unsigned int hit = 0x80000000u;
    optixTrace(params.traversable,
               origin,
               direction,
               tMin,
               tMax,
               0.0f,
               OptixVisibilityMask(255),
               OPTIX_RAY_FLAG_DISABLE_CLOSESTHIT | OPTIX_RAY_FLAG_ENFORCE_ANYHIT |
                   OPTIX_RAY_FLAG_TERMINATE_ON_FIRST_HIT,
               0,
               1,
               0,
               hit);
    return hit;
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
    const unsigned int payload = optixGetPayload_0();
    if (payload == 0x80000000u)
    {
        optixSetPayload_0(0u);
        return;
    }

    const float3 direction = optixGetWorldRayDirection();
    const float3 color = SkyColor(direction);
    optixSetPayload_0(PackColor(color.x, color.y, color.z));
}

extern "C" __global__ void __anyhit__primary()
{
    const unsigned int payload = optixGetPayload_0();
    if (payload == 0x80000000u)
    {
        optixSetPayload_0(1u);
        optixTerminateRay();
    }
}

extern "C" __global__ void __closesthit__primary()
{
    if (params.integrator == 1)
    {
        const int spp = max(params.spp, 1);
        const uint3 launchIndex = optixGetLaunchIndex();
        unsigned int rngState = Hash32((launchIndex.x + 1u) * 73856093u ^ (launchIndex.y + 1u) * 19349663u);

        const float3 rayOrigin = optixGetWorldRayOrigin();
        const float3 rayDirection = Normalize3(optixGetWorldRayDirection());
        const float tHit = optixGetRayTmax();
        const float3 hitPoint =
            make_float3(rayOrigin.x + rayDirection.x * tHit,
                        rayOrigin.y + rayDirection.y * tHit,
                        rayOrigin.z + rayDirection.z * tHit);
        float3 normal = Normalize3(make_float3(-rayDirection.x, -rayDirection.y, -rayDirection.z));
        const unsigned int hitKind = optixGetHitKind();
        if (hitKind == OPTIX_HIT_KIND_TRIANGLE_FRONT_FACE ||
            hitKind == OPTIX_HIT_KIND_TRIANGLE_BACK_FACE)
        {
            const unsigned int instanceId = optixGetInstanceId();
            if (params.instanceGeomRefs != 0ull && instanceId < (unsigned int)params.instanceGeomRefCount)
            {
                const LaunchParams::InstanceGeomRef *refs =
                    reinterpret_cast<const LaunchParams::InstanceGeomRef *>(params.instanceGeomRefs);
                const LaunchParams::InstanceGeomRef ref = refs[instanceId];
                const int primitiveIndex = int(optixGetPrimitiveIndex());
                const int indexBase = primitiveIndex * 3;
                const float3 *positions = reinterpret_cast<const float3 *>(ref.positions);
                const int *indices = reinterpret_cast<const int *>(ref.indices);
                if (positions != nullptr && indices != nullptr && indexBase + 2 < ref.numIndices)
                {
                    const int i0 = indices[indexBase + 0];
                    const int i1 = indices[indexBase + 1];
                    const int i2 = indices[indexBase + 2];
                    if (i0 >= 0 && i0 < ref.numPositions && i1 >= 0 && i1 < ref.numPositions &&
                        i2 >= 0 && i2 < ref.numPositions)
                    {
                        const float3 p0 = positions[i0];
                        const float3 p1 = positions[i1];
                        const float3 p2 = positions[i2];
                        const float3 edge01 = make_float3(p1.x - p0.x, p1.y - p0.y, p1.z - p0.z);
                        const float3 edge02 = make_float3(p2.x - p0.x, p2.y - p0.y, p2.z - p0.z);
                        const float3 localNormal = Normalize3(Cross3(edge01, edge02));
                        const float3 worldNormal =
                            Normalize3(optixTransformNormalFromObjectToWorldSpace(localNormal));
                        normal = FaceForward(worldNormal, rayDirection);
                    }
                }
            }
        }
        float3 tangent, bitangent;
        BuildOrthonormalBasis(normal, tangent, bitangent);

        int visible = 0;
        for (int i = 0; i < spp; i++)
        {
            const float u1 = Random01(rngState);
            const float u2 = Random01(rngState);
            const float3 local = SampleCosineHemisphere(u1, u2);
            const float3 sampleDir = Normalize3(make_float3(tangent.x * local.x + bitangent.x * local.y + normal.x * local.z,
                                                            tangent.y * local.x + bitangent.y * local.y + normal.y * local.z,
                                                            tangent.z * local.x + bitangent.z * local.y + normal.z * local.z));
            const float3 sampleOrigin = make_float3(hitPoint.x + normal.x * params.aoBias,
                                                    hitPoint.y + normal.y * params.aoBias,
                                                    hitPoint.z + normal.z * params.aoBias);
            const unsigned int occluded =
                TraceOcclusion(sampleOrigin, sampleDir, params.aoBias, params.aoMaxDistance);
            if (occluded == 0u)
            {
                visible++;
            }
        }

        const float ao = float(visible) / float(spp);
        const float3 aoColor = make_float3(ao, ao, ao);
        optixSetPayload_0(PackColor(aoColor.x, aoColor.y, aoColor.z));
        return;
    }

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
