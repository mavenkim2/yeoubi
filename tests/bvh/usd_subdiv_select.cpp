#include "usd_subdiv_select.h"

#include <pxr/usd/usd/primRange.h>
#include <pxr/usd/usdGeom/mesh.h>
#include <pxr/usd/usdGeom/tokens.h>

bool SelectLargestCatmullClarkMesh(const pxr::UsdStageRefPtr &stage, SelectedSubdivMesh &meshOut)
{
    if (!stage)
    {
        return false;
    }

    const pxr::Usd_PrimFlagsConjunction flags =
        pxr::UsdPrimIsActive && pxr::UsdPrimIsLoaded && !pxr::UsdPrimIsAbstract;
    const pxr::Usd_PrimFlagsPredicate predicate(flags);

    bool found = false;
    size_t maxFaceCount = 0;

    for (pxr::UsdPrim prim : stage->Traverse(predicate))
    {
        if (!prim.IsA<pxr::UsdGeomMesh>())
        {
            continue;
        }

        pxr::UsdGeomMesh mesh(prim);
        pxr::TfToken scheme;
        if (!mesh.GetSubdivisionSchemeAttr().Get(&scheme))
        {
            continue;
        }
        if (scheme != pxr::UsdGeomTokens->catmullClark)
        {
            continue;
        }

        pxr::VtVec3fArray points;
        pxr::VtIntArray faceCounts;
        pxr::VtIntArray faceIndices;
        if (!mesh.GetPointsAttr().Get(&points) || !mesh.GetFaceVertexCountsAttr().Get(&faceCounts) ||
            !mesh.GetFaceVertexIndicesAttr().Get(&faceIndices))
        {
            continue;
        }

        const size_t faceCount = faceCounts.size();
        if (!found || faceCount > maxFaceCount)
        {
            meshOut.path = prim.GetPath();
            meshOut.points = std::move(points);
            meshOut.faceVertexCounts = std::move(faceCounts);
            meshOut.faceVertexIndices = std::move(faceIndices);
            meshOut.subdivisionScheme = scheme.GetString();
            maxFaceCount = faceCount;
            found = true;
        }
    }

    return found;
}
