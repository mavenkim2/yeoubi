#pragma once

#include "kernel/state/types.h"
#include <optix_device.h>

namespace ybi
{

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

} // namespace ybi
