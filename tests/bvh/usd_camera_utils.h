#pragma once

#include <pxr/base/gf/vec3f.h>
#include <pxr/usd/sdf/path.h>
#include <pxr/usd/usd/primRange.h>
#include <pxr/usd/usd/stage.h>
#include <pxr/usd/usdGeom/camera.h>
#include <pxr/usd/usdGeom/xformCache.h>

#include <algorithm>
#include <limits>

struct UsdCameraInfo
{
    bool found = false;
    pxr::SdfPath path;
    pxr::GfVec3f worldPosition;
    pxr::GfVec3f meshCenter;
    float distanceToMeshCenter = 0.0f;
};

static inline pxr::GfVec3f ComputeMeshCenter(const pxr::VtVec3fArray &points)
{
    if (points.empty())
    {
        return pxr::GfVec3f(0.0f);
    }

    pxr::GfVec3f bbMin(std::numeric_limits<float>::max());
    pxr::GfVec3f bbMax(std::numeric_limits<float>::lowest());
    for (const pxr::GfVec3f &p : points)
    {
        bbMin[0] = std::min(bbMin[0], p[0]);
        bbMin[1] = std::min(bbMin[1], p[1]);
        bbMin[2] = std::min(bbMin[2], p[2]);
        bbMax[0] = std::max(bbMax[0], p[0]);
        bbMax[1] = std::max(bbMax[1], p[1]);
        bbMax[2] = std::max(bbMax[2], p[2]);
    }
    return (bbMin + bbMax) * 0.5f;
}

static inline bool GetClosestUsdCameraInfo(const pxr::UsdStageRefPtr &stage,
                                           const pxr::VtVec3fArray &meshPoints,
                                           UsdCameraInfo &cameraInfoOut)
{
    cameraInfoOut = {};
    if (!stage)
    {
        return false;
    }

    const pxr::GfVec3f center = ComputeMeshCenter(meshPoints);
    const pxr::Usd_PrimFlagsConjunction flags =
        pxr::UsdPrimIsActive && pxr::UsdPrimIsLoaded && !pxr::UsdPrimIsAbstract;
    const pxr::Usd_PrimFlagsPredicate predicate(flags);

    pxr::UsdGeomXformCache xformCache(pxr::UsdTimeCode::Default());
    float closestDistance = std::numeric_limits<float>::max();
    for (const pxr::UsdPrim &prim : stage->Traverse(predicate))
    {
        if (!prim.IsA<pxr::UsdGeomCamera>())
        {
            continue;
        }

        const pxr::GfMatrix4d localToWorld = xformCache.GetLocalToWorldTransform(prim);
        const pxr::GfVec3d positionD = localToWorld.Transform(pxr::GfVec3d(0.0, 0.0, 0.0));
        const pxr::GfVec3f cameraPosition{
            float(positionD[0]), float(positionD[1]), float(positionD[2])};
        const float distanceToCenter = (cameraPosition - center).GetLength();

        if (distanceToCenter < closestDistance)
        {
            closestDistance = distanceToCenter;
            cameraInfoOut.found = true;
            cameraInfoOut.path = prim.GetPath();
            cameraInfoOut.worldPosition = cameraPosition;
            cameraInfoOut.meshCenter = center;
            cameraInfoOut.distanceToMeshCenter = distanceToCenter;
        }
    }

    return cameraInfoOut.found;
}
