#pragma once

#include "core/vec3f.h"
#include "device/cuda_compat.h"
#include "kernel/state/types.h"

#if defined(__CUDACC__) && defined(OPTIX_VERSION)
#include "kernel/backend/optix/ray_tracing.h"
#elif defined(__CUDACC__)
#error "Ray tracing requires OptiX; CUDA-only ray tracing is not implemented."
#else
#error "Ray tracing requires OptiX device compilation; no supported ray-tracing backend is active."
#endif

namespace ybi
{

YBI_GPU IntersectionInfo TraceRay(uint64_t bvhHandle,
                                  const Vec3f &rayOrigin,
                                  const Vec3f &rayDirection,
                                  float tMin,
                                  float tMax)
{
    return TraceOptixRay(bvhHandle, rayOrigin, rayDirection, tMin, tMax);
}

YBI_GPU bool TraceVisibility(uint64_t bvhHandle,
                             const Vec3f &rayOrigin,
                             const Vec3f &rayDirection,
                             float tMax)
{
    return TraceOptixVisibility(bvhHandle, rayOrigin, rayDirection, tMax);
}

} // namespace ybi
