#include "usd_subdiv_select.h"

#include <pxr/usd/usd/primRange.h>
#include <pxr/usd/usdGeom/mesh.h>
#include <pxr/usd/usdGeom/tokens.h>

static pxr::TfToken GetTokenAttrOrDefault(const pxr::UsdAttribute &attr, const pxr::TfToken &fallback)
{
    pxr::TfToken value = fallback;
    if (attr)
    {
        attr.Get(&value);
    }
    return value;
}

static pxr::TfToken GetCreaseMethodOrDefault(const pxr::UsdPrim &prim)
{
    pxr::TfToken value = pxr::TfToken("uniform");
    const pxr::UsdAttribute attr = prim.GetAttribute(pxr::TfToken("creaseMethod"));
    if (attr && attr.Get(&value))
    {
        return value;
    }
    const pxr::UsdAttribute namespaced = prim.GetAttribute(pxr::TfToken("osd:creaseMethod"));
    if (namespaced && namespaced.Get(&value))
    {
        return value;
    }
    return value;
}

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
        pxr::TfToken vtxBoundary;
        pxr::TfToken fvarLinear;
        pxr::TfToken triSub;
        if (!mesh.GetSubdivisionSchemeAttr().Get(&scheme))
        {
            continue;
        }
        if (scheme != pxr::UsdGeomTokens->catmullClark)
        {
            continue;
        }
        vtxBoundary = GetTokenAttrOrDefault(mesh.GetInterpolateBoundaryAttr(),
                                            pxr::UsdGeomTokens->edgeAndCorner);
        fvarLinear = GetTokenAttrOrDefault(mesh.GetFaceVaryingLinearInterpolationAttr(),
                                           pxr::UsdGeomTokens->cornersPlus1);
        triSub = GetTokenAttrOrDefault(mesh.GetTriangleSubdivisionRuleAttr(),
                                       pxr::UsdGeomTokens->catmullClark);
        const pxr::TfToken creasing = GetCreaseMethodOrDefault(prim);

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
            meshOut.vertexBoundaryInterpolation = vtxBoundary.GetString();
            meshOut.fvarLinearInterpolation = fvarLinear.GetString();
            meshOut.creasingMethod = creasing.GetString();
            meshOut.triangleSubdivision = triSub.GetString();
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
