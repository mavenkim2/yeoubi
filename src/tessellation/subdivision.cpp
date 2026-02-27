#include "tessellation/subdivision.h"

#include "bvh/usd_subdiv_select.h"
#include "tessellation/edge_map_validation.h"
#include "tessellation/subdivision_patch_types.h"
#include "util/assert.h"

#include <opensubdiv/far/patchMap.h>
#include <opensubdiv/far/patchTableFactory.h>
#include <opensubdiv/far/primvarRefiner.h>
#include <opensubdiv/far/ptexIndices.h>
#include <opensubdiv/far/topologyDescriptor.h>
#include <opensubdiv/far/topologyRefinerFactory.h>
#include <pxr/base/gf/vec2f.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <string>
#include <unordered_map>
#include <vector>

using namespace OpenSubdiv;

YBI_NAMESPACE_BEGIN

struct LimitEvalVertex
{
    pxr::GfVec3f p = pxr::GfVec3f(0.0f);

    void Clear()
    {
        p = pxr::GfVec3f(0.0f);
    }

    void AddWithWeight(const LimitEvalVertex &src, float w)
    {
        p += src.p * w;
    }
};

// clang-format off
#include "tessellation/tessellation_adaptive_limit_eval.h"
#include "tessellation/tessellation_adaptive_edge_ops.h"
#include "tessellation/tessellation_adaptive_patch_build.h"
#include "tessellation/tessellation_adaptive_util.h"
#include "tessellation/tessellation_adaptive_patch_quad_write.h"
#include "tessellation/tessellation_adaptive_obj_write.h"
#include "tessellation/tessellation_adaptive_mesh_build.h"
// clang-format on

static SelectedSubdivMesh ConvertToSelectedSubdivMesh(const SubdivisionMesh &mesh,
                                                      const SubdivisionRunOptions &options)
{
    SelectedSubdivMesh out = {};
    out.path = pxr::SdfPath::EmptyPath();

    out.points.resize(mesh.vertices.size());
    for (size_t i = 0; i < mesh.vertices.size(); ++i)
    {
        const float3 p = mesh.vertices[i];
        out.points[i] = pxr::GfVec3f(p.x, p.y, p.z);
    }

    out.faceVertexCounts.resize(mesh.vertsPerFace.size());
    for (size_t i = 0; i < mesh.vertsPerFace.size(); ++i)
    {
        out.faceVertexCounts[i] = mesh.vertsPerFace[i];
    }

    out.faceVertexIndices.resize(mesh.indices.size());
    for (size_t i = 0; i < mesh.indices.size(); ++i)
    {
        out.faceVertexIndices[i] = mesh.indices[i];
    }

    out.cornerIndices.resize(mesh.cornerIndices.size());
    for (size_t i = 0; i < mesh.cornerIndices.size(); ++i)
    {
        out.cornerIndices[i] = mesh.cornerIndices[i];
    }

    out.cornerSharpnesses.resize(mesh.cornerSharpnesses.size());
    for (size_t i = 0; i < mesh.cornerSharpnesses.size(); ++i)
    {
        out.cornerSharpnesses[i] = mesh.cornerSharpnesses[i];
    }

    out.creaseIndices.resize(mesh.creaseIndices.size());
    for (size_t i = 0; i < mesh.creaseIndices.size(); ++i)
    {
        out.creaseIndices[i] = mesh.creaseIndices[i];
    }

    out.creaseLengths.resize(mesh.creaseLengths.size());
    for (size_t i = 0; i < mesh.creaseLengths.size(); ++i)
    {
        out.creaseLengths[i] = mesh.creaseLengths[i];
    }

    out.creaseSharpnesses.resize(mesh.creaseSharpnesses.size());
    for (size_t i = 0; i < mesh.creaseSharpnesses.size(); ++i)
    {
        out.creaseSharpnesses[i] = mesh.creaseSharpnesses[i];
    }

    out.holeIndices.resize(mesh.holeIndices.size());
    for (size_t i = 0; i < mesh.holeIndices.size(); ++i)
    {
        out.holeIndices[i] = mesh.holeIndices[i];
    }

    out.subdivisionScheme = options.subdivisionScheme;
    switch (mesh.interpolationRule)
    {
        case BOUNDARY_INTERPOLATION_NONE: out.vertexBoundaryInterpolation = "none"; break;
        case BOUNDARY_INTERPOLATION_EDGE: out.vertexBoundaryInterpolation = "edgeOnly"; break;
        case BOUNDARY_INTERPOLATION_EDGE_AND_CORNER:
        default: out.vertexBoundaryInterpolation = "edgeAndCorner"; break;
    }
    switch (mesh.fvarLinearInterpolation)
    {
        case FVAR_LINEAR_NONE: out.fvarLinearInterpolation = "none"; break;
        case FVAR_LINEAR_CORNERS_ONLY: out.fvarLinearInterpolation = "cornersOnly"; break;
        case FVAR_LINEAR_CORNERS_PLUS1: out.fvarLinearInterpolation = "cornersPlus1"; break;
        case FVAR_LINEAR_CORNERS_PLUS2: out.fvarLinearInterpolation = "cornersPlus2"; break;
        case FVAR_LINEAR_BOUNDARIES: out.fvarLinearInterpolation = "boundaries"; break;
        case FVAR_LINEAR_ALL: out.fvarLinearInterpolation = "all"; break;
        default: out.fvarLinearInterpolation = "cornersPlus1"; break;
    }
    out.creasingMethod = options.creasingMethod;
    out.triangleSubdivision = options.triangleSubdivision;

    return out;
}

static pxr::GfVec3f ComputeMeshCenter(const pxr::VtVec3fArray &points)
{
    if (points.empty())
    {
        return pxr::GfVec3f(0.0f, 0.0f, 0.0f);
    }

    pxr::GfVec3f sum(0.0f, 0.0f, 0.0f);
    for (const auto &p : points)
    {
        sum += p;
    }
    return sum * (1.0f / float(points.size()));
}

bool SubdivideAdaptive(const SubdivisionMesh &mesh,
                      const SubdivisionRunOptions &options,
                      SubdivisionRunResult *outResult)
{
    if (!outResult)
    {
        return false;
    }

    SelectedSubdivMesh m = ConvertToSelectedSubdivMesh(mesh, options);

    CreasePairs creases = BuildCreasePairs(m);
    Far::TopologyDescriptor d = {};
    d.numVertices = int(m.points.size());
    d.numFaces = int(m.faceVertexCounts.size());
    d.numVertsPerFace = m.faceVertexCounts.data();
    d.vertIndicesPerFace = m.faceVertexIndices.data();
    d.numCreases = int(creases.weights.size());
    d.creaseVertexIndexPairs = creases.pairs.data();
    d.creaseWeights = creases.weights.data();
    d.numCorners = int(m.cornerIndices.size());
    d.cornerVertexIndices = m.cornerIndices.data();
    d.cornerWeights = m.cornerSharpnesses.data();
    d.numHoles = int(m.holeIndices.size());
    d.holeIndices = m.holeIndices.data();

    Sdc::Options sdcOptions;
    sdcOptions.SetVtxBoundaryInterpolation(VtxBoundaryFromString(m.vertexBoundaryInterpolation));
    sdcOptions.SetFVarLinearInterpolation(FVarLinearFromString(m.fvarLinearInterpolation));
    sdcOptions.SetCreasingMethod(CreasingMethodFromString(m.creasingMethod));
    sdcOptions.SetTriangleSubdivision(TriangleSubFromString(m.triangleSubdivision));

    Far::TopologyRefinerFactory<Far::TopologyDescriptor>::Options o(
        SchemeFromString(m.subdivisionScheme), sdcOptions);

    SubdivisionEdgeMap edgeMap;
    edgeMap.reserve(m.faceVertexIndices.size());
    const int edgesWithOver2Faces = CountNonManifoldEdges(edgeMap);

    Far::TopologyRefiner *refiner =
        Far::TopologyRefinerFactory<Far::TopologyDescriptor>::Create(d, o);
    if (!refiner)
    {
        return false;
    }

    int nextGeneratedVertexId = int(m.points.size());
    const std::vector<SubdivisionPatch> patches =
        BuildSubdivisionPatches(m, *refiner, edgeMap, nextGeneratedVertexId);

    const EdgeMapChecks edgeChecks = RunEdgeMapChecks(m, patches, edgeMap);
    if (!edgeChecks.ok)
    {
        delete refiner;
        return false;
    }

    Far::TopologyRefiner::AdaptiveOptions adaptiveOptions(options.level);
    refiner->RefineAdaptive(adaptiveOptions);
    if (refiner->GetMaxLevel() < options.level)
    {
        delete refiner;
        return false;
    }

    Far::PatchTableFactory::Options patchOptions(options.level);
    patchOptions.endCapType = Far::PatchTableFactory::Options::ENDCAP_GREGORY_BASIS;
    patchOptions.useInfSharpPatch = (d.numCreases > 0) || (d.numCorners > 0);
    patchOptions.generateFVarLegacyLinearPatches = false;
    const Far::PatchTable *patchTable = Far::PatchTableFactory::Create(*refiner, patchOptions);
    if (!patchTable)
    {
        delete refiner;
        return false;
    }

    Far::PatchMap patchMap(*patchTable);
    const std::vector<LimitEvalVertex> limitValues =
        BuildLimitEvalVertices(*refiner, *patchTable, m.points);

    const pxr::GfVec3f eye = options.eye;
    const pxr::GfVec3f lookAt = options.lookAt;

    int tmaxComputedEdges = 0;
    const std::vector<SubdivisionPatch> splitPatches =
        DiagSplitPatches(m,
                         patches,
                         edgeMap,
                         nextGeneratedVertexId,
                         patchMap,
                         *patchTable,
                         limitValues,
                         options.sampleSteps,
                         options.pixelSpacing,
                         options.splitThreshold,
                         eye,
                         lookAt,
                         options.viewportWidth,
                         options.viewportHeight,
                         options.verticalFovDegrees,
                         &tmaxComputedEdges);

    int missingStoredPatchParams = 0;
    int badUniformFactor = 0;
    if (!VerifyInitializedUniformEdgesHaveStoredPatchParams(
            edgeMap, &missingStoredPatchParams, &badUniformFactor))
    {
        delete patchTable;
        delete refiner;
        return false;
    }

    int patchQuadVerts = 0;
    int patchQuadCount = 0;
    if (!options.patchQuadObjPath.empty())
    {
        if (!WriteLeafPatchCornerQuadsObj(splitPatches,
                                          edgeMap,
                                          patchMap,
                                          *patchTable,
                                          limitValues,
                                          options.patchQuadObjPath,
                                          &patchQuadVerts,
                                          &patchQuadCount))
        {
            delete patchTable;
            delete refiner;
            return false;
        }
    }

    int innerGridVerts = 0;
    int innerGridTris = 0;
    if (!BuildLeafPatchStitchedMesh(splitPatches,
                                    edgeMap,
                                    patchMap,
                                    *patchTable,
                                    limitValues,
                                    nextGeneratedVertexId,
                                    &outResult->mesh,
                                    &outResult->trianglePatchFaceIds,
                                    &outResult->triangleCoarseFaceIds,
                                    &outResult->trianglePtexFaceIds,
                                    &outResult->triangleQuadrants,
                                    &innerGridVerts,
                                    &innerGridTris))
    {
        delete patchTable;
        delete refiner;
        return false;
    }

    outResult->refinedMaxLevel = refiner->GetMaxLevel();
    outResult->totalPatches = patchTable->GetNumPatchesTotal();
    outResult->subdivisionPatchCount = patches.size();
    outResult->diagSplitPatchCount = splitPatches.size();
    outResult->generatedVertexCount = nextGeneratedVertexId - int(m.points.size());
    outResult->midpointEdges = CountEdgesWithMidpointVertex(edgeMap);
    outResult->edgeTMaxComputed = tmaxComputedEdges;
    outResult->totalComputedEdges = CountEdgesWithComputedTMax(edgeMap);
    outResult->patchQuadVerts = patchQuadVerts;
    outResult->patchQuadCount = patchQuadCount;
    outResult->controlCageUniqueEdges = edgeMap.size();
    outResult->boundaryEdges = CountBoundaryEdges(edgeMap);
    outResult->controlCageEdgesWithOver2Faces = edgesWithOver2Faces;

    delete patchTable;
    delete refiner;
    return true;
}

YBI_NAMESPACE_END
