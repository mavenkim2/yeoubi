#pragma once

#include "device/cuda_compat.h"
#include "kernel/geometry/triangle_info.h"
#include "kernel/shading/shading.h"
#include "kernel/state/types.h"

namespace ybi
{

// Improved Shader and Texture Level of Detail Using Ray Cones (
YBI_GPU void PropagateRayCone(RayCone &cone, float surfaceSpreadAngle, float tHit)
{
    cone.width += cone.spreadAngle * tHit;
    cone.spreadAngle += surfaceSpreadAngle;
}

YBI_GPU Vec3f CalculateEdgeCurvatures(const Vec3f &p0,
                                      const Vec3f &p1,
                                      const Vec3f &p2,
                                      const Vec3f &n0,
                                      const Vec3f &n1,
                                      const Vec3f &n2)
{
    Vec3f p1p0 = p1 - p0;
    Vec3f p2p1 = p2 - p1;
    Vec3f p0p2 = p0 - p2;
    float k01 = Dot(n1 - n0, p1p0) / Dot(p1p0, p1p0);
    float k12 = Dot(n2 - n1, p2p1) / Dot(p2p1, p2p1);
    float k20 = Dot(n0 - n2, p0p2) / Dot(p0p2, p0p2);
    return Vec3f(k01, k12, k20);
}

YBI_GPU float CalculateSurfaceSpreadAngle(const RayCone &cone, float k, Vec3f d, Vec3f n)
{
    float bc = -2.f * k * cone.width / Dot(d, n);
    return bc;
}

YBI_GPU float CalculatePrimaryHitUnifiedSurfaceSpreadAngle(
    RayCone &cone, Vec3f d, Vec3f n, Vec3f p0, Vec3f p1, Vec3f p2, Vec3f n0, Vec3f n1, Vec3f n2)
{
    Vec3f k = CalculateEdgeCurvatures(p0, p1, p2, n0, n1, n2);
    float minK = Min(k.x, Min(k.y, k.z));
    float maxK = Max(k.x, Max(k.y, k.z));
    float spreadMinK = CalculateSurfaceSpreadAngle(cone, minK, d, n);
    float spreadMaxK = CalculateSurfaceSpreadAngle(cone, maxK, d, n);
    float test0 = Abs(cone.spreadAngle + spreadMinK);
    float test1 = Abs(cone.spreadAngle + spreadMaxK);
    float result = test0 >= test1 ? spreadMinK : spreadMaxK;
    return result;
}

YBI_GPU void UpdatePrimaryRayCone(RayCone &cone,
                                  const Vec3f &rayDirection,
                                  float rayT,
                                  const TriangleInfo &triangleInfo,
                                  const TriangleSurfaceInfo &surfaceInfo)
{
    const Vec3f fallbackNormal = surfaceInfo.geometricNormal;
    const Vec3f n0 = triangleInfo.hasNormals ? triangleInfo.normals[0] : fallbackNormal;
    const Vec3f n1 = triangleInfo.hasNormals ? triangleInfo.normals[1] : fallbackNormal;
    const Vec3f n2 = triangleInfo.hasNormals ? triangleInfo.normals[2] : fallbackNormal;
    const float surfaceSpreadAngle =
        CalculatePrimaryHitUnifiedSurfaceSpreadAngle(cone,
                                                     rayDirection,
                                                     surfaceInfo.shadingNormal,
                                                     triangleInfo.positions[0],
                                                     triangleInfo.positions[1],
                                                     triangleInfo.positions[2],
                                                     n0,
                                                     n1,
                                                     n2);
    PropagateRayCone(cone, surfaceSpreadAngle, rayT);
}

YBI_GPU float CalculateSecondaryHitSurfaceSpreadAngle(const Vec3f &d,
                                                      const Vec3f &n,
                                                      const Vec3f &p0,
                                                      const Vec3f &p1,
                                                      const Vec3f &p2,
                                                      const Vec3f &n0,
                                                      const Vec3f &n1,
                                                      const Vec3f &n2)
{
    Vec3f k = CalculateEdgeCurvatures(p0, p1, p2, n0, n1, n2);
    float result = (k.x + k.y + k.z) / 3.f;
    return result;
}

YBI_GPU void UpdateSecondaryRayCone(RayCone &cone,
                                    const Vec3f &rayDirection,
                                    float rayT,
                                    const TriangleInfo &triangleInfo,
                                    const TriangleSurfaceInfo &surfaceInfo)
{
    float surfaceSpreadAngle = CalculateSecondaryHitSurfaceSpreadAngle(rayDirection,
                                                                       surfaceInfo.shadingNormal,
                                                                       triangleInfo.positions[0],
                                                                       triangleInfo.positions[1],
                                                                       triangleInfo.positions[2],
                                                                       triangleInfo.normals[0],
                                                                       triangleInfo.normals[1],
                                                                       triangleInfo.normals[2]);
    PropagateRayCone(cone, surfaceSpreadAngle, rayT);
}

YBI_GPU void UpdateRayCone(RayCone &cone,
                           const Vec3f &rayDirection,
                           const float rayT,
                           const ShadingInfo &info,
                           const uint32_t depth)
{
    // Primary hit
    float surfaceSpreadAngle;
    if (depth == 0)
    {
        surfaceSpreadAngle = CalculatePrimaryHitUnifiedSurfaceSpreadAngle(cone,
                                                                          rayDirection,
                                                                          info.shadingNormal,
                                                                          info.positions[0],
                                                                          info.positions[1],
                                                                          info.positions[2],
                                                                          info.normals[0],
                                                                          info.normals[1],
                                                                          info.normals[2]);
    }
    else
    {
        surfaceSpreadAngle = CalculateSecondaryHitSurfaceSpreadAngle(rayDirection,
                                                                     info.shadingNormal,
                                                                     info.positions[0],
                                                                     info.positions[1],
                                                                     info.positions[2],
                                                                     info.normals[0],
                                                                     info.normals[1],
                                                                     info.normals[2]);
    }
    PropagateRayCone(cone, surfaceSpreadAngle, rayT);
}

} // namespace ybi
