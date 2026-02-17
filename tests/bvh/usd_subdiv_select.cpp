#include "usd_subdiv_select.h"

#include <pxr/usd/usd/primRange.h>
#include <pxr/usd/usdGeom/mesh.h>
#include <pxr/usd/usdGeom/tokens.h>

static bool HasValidCreaseData(const pxr::VtIntArray &creaseLengths,
                               const pxr::VtFloatArray &creaseSharpnesses)
{
    const size_t creaseCount = std::min(creaseLengths.size(), creaseSharpnesses.size());
    for (size_t i = 0; i < creaseCount; i++)
    {
        if (creaseLengths[i] >= 2 && creaseSharpnesses[i] > 0.0f)
        {
            return true;
        }
    }
    return false;
}

static bool SelectLargestCatmullClarkMeshImpl(const pxr::UsdStageRefPtr &stage,
                                              SelectedSubdivMesh &meshOut,
                                              bool requireCreases)
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
        pxr::VtIntArray cornerIndices;
        pxr::VtFloatArray cornerSharpnesses;
        pxr::VtIntArray creaseIndices;
        pxr::VtIntArray creaseLengths;
        pxr::VtFloatArray creaseSharpnesses;
        pxr::VtIntArray holeIndices;
        if (!mesh.GetPointsAttr().Get(&points) || !mesh.GetFaceVertexCountsAttr().Get(&faceCounts) ||
            !mesh.GetFaceVertexIndicesAttr().Get(&faceIndices))
        {
            continue;
        }
        mesh.GetCornerIndicesAttr().Get(&cornerIndices);
        mesh.GetCornerSharpnessesAttr().Get(&cornerSharpnesses);
        mesh.GetCreaseIndicesAttr().Get(&creaseIndices);
        mesh.GetCreaseLengthsAttr().Get(&creaseLengths);
        mesh.GetCreaseSharpnessesAttr().Get(&creaseSharpnesses);
        mesh.GetHoleIndicesAttr().Get(&holeIndices);
        if (requireCreases && !HasValidCreaseData(creaseLengths, creaseSharpnesses))
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
            meshOut.cornerIndices = std::move(cornerIndices);
            meshOut.cornerSharpnesses = std::move(cornerSharpnesses);
            meshOut.creaseIndices = std::move(creaseIndices);
            meshOut.creaseLengths = std::move(creaseLengths);
            meshOut.creaseSharpnesses = std::move(creaseSharpnesses);
            meshOut.holeIndices = std::move(holeIndices);
            meshOut.subdivisionScheme = scheme.GetString();
            maxFaceCount = faceCount;
            found = true;
        }
    }

    return found;
}

bool SelectLargestCatmullClarkMesh(const pxr::UsdStageRefPtr &stage, SelectedSubdivMesh &meshOut)
{
    return SelectLargestCatmullClarkMeshImpl(stage, meshOut, false);
}

bool SelectLargestCatmullClarkMeshWithCreases(const pxr::UsdStageRefPtr &stage,
                                              SelectedSubdivMesh &meshOut)
{
    return SelectLargestCatmullClarkMeshImpl(stage, meshOut, true);
}
