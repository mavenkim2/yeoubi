#include "kernel/integrator/direct_lighting.h"
#include "kernel/integrator/indirect_lighting.h"
#include "kernel/integrator/intersect.h"
#include "kernel/integrator/shade_hit.h"
#include "kernel/launch/wavefront_params.h"
#include <optix_device.h>

namespace ybi
{

extern "C"
{
    __constant__ WavefrontParams params;
}

YBI_GPU IntersectionInfo *GetOptixIntersectionPayload()
{
    const uint64_t ptr = uint64_t(optixGetPayload_0()) |
                         (uint64_t(optixGetPayload_1()) << 32u);
    return reinterpret_cast<IntersectionInfo *>(ptr);
}

extern "C" __global__ void __miss__kernel()
{
    GetOptixIntersectionPayload()->type = INTERSECTION_TYPE_NONE;
}

extern "C" __global__ void __closesthit__kernel()
{
    IntersectionInfo *isect = GetOptixIntersectionPayload();
    if (optixIsTriangleHit())
    {
        const float2 bary = optixGetTriangleBarycentrics();
        float objectToWorld[12];
        optixGetObjectToWorldTransformMatrix(objectToWorld);

        isect->t = optixGetRayTmax();
        isect->id = optixGetInstanceId();
        isect->primIndex = optixGetPrimitiveIndex();
        isect->clusterId = ~0u;
#if OPTIX_VERSION >= 90000
        isect->clusterId = optixGetClusterId();
#endif
        isect->bary = Vec2f(bary.x, bary.y);
        isect->objectToWorld = AffineTransform(objectToWorld);
        isect->type = INTERSECTION_TYPE_TRIANGLE;
    }
}

extern "C" __global__ void __raygen__intersect()
{
    const uint32_t threadIndex = optixGetLaunchIndex().x;
    if (threadIndex >= params.numSampleIndices)
        return;

    // Initial camera batches occupy contiguous state slots.
    const uint32_t sampleIndex =
        params.sampleIndices ? params.sampleIndices[threadIndex] : threadIndex;
    Intersect(params, &params.sampleState, sampleIndex);
}

extern "C" __global__ void __raygen__shade_hit()
{
    const uint32_t threadIndex = optixGetLaunchIndex().x;
    if (threadIndex >= params.numSampleIndices)
        return;

    const uint32_t sampleIndex = params.sampleIndices[threadIndex];
    ShadeHit(params, &params.sampleState, sampleIndex);
}

extern "C" __global__ void __raygen__direct_lighting()
{
    const uint32_t threadIndex = optixGetLaunchIndex().x;
    if (threadIndex >= params.numSampleIndices)
        return;

    const uint32_t sampleIndex = params.sampleIndices[threadIndex];
    DirectLighting(params, &params.sampleState, sampleIndex);
}

extern "C" __global__ void __raygen__indirect_lighting()
{
    const uint32_t threadIndex = optixGetLaunchIndex().x;
    if (threadIndex >= params.numSampleIndices)
        return;

    const uint32_t sampleIndex = params.sampleIndices[threadIndex];
    IndirectLighting(params, &params.sampleState, sampleIndex);
}

} // namespace ybi
