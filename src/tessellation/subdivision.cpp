#include "tessellation/subdivision.h"

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
#include <cstdio>
#include <limits>
#include <string>
#include <unordered_map>
#include <vector>

using namespace OpenSubdiv;

YBI_NAMESPACE_BEGIN

struct LimitEvalVertex
{
    float3 p = make_float3(0.0f);

    void Clear()
    {
        p = make_float3(0.0f);
    }

    void AddWithWeight(const LimitEvalVertex &src, float w)
    {
        p += src.p * w;
    }
};

struct LimitEvalFVar2
{
    float2 uv = make_float2(0.0f);

    void Clear()
    {
        uv = make_float2(0.0f);
    }

    void AddWithWeight(const LimitEvalFVar2 &src, float w)
    {
        uv += src.uv * w;
    }
};

// clang-format off
#include "tessellation/tessellation_adaptive_limit_eval.h"
#include "tessellation/tessellation_adaptive_edge_ops.h"
#include "tessellation/tessellation_adaptive_patch_build.h"
#include "tessellation/tessellation_adaptive_util.h"
#include "tessellation/tessellation_adaptive_patch_quad_write.h"
#include "tessellation/tessellation_adaptive_mesh_build.h"
// clang-format on

bool SubdivideAdaptive(const SubdivisionMesh &mesh,
                      const SubdivisionRunOptions &options,
                      SubdivisionRunResult *outResult)
{
    if (!outResult)
    {
        return false;
    }

    CreasePairs creases = BuildCreasePairs(mesh);
    Far::TopologyDescriptor d = {};
    d.numVertices = int(mesh.vertices.size());
    d.numFaces = int(mesh.vertsPerFace.size());
    d.numVertsPerFace = mesh.vertsPerFace.data();
    d.vertIndicesPerFace = mesh.indices.data();
    d.numCreases = int(creases.weights.size());
    d.creaseVertexIndexPairs = creases.pairs.data();
    d.creaseWeights = creases.weights.data();
    d.numCorners = int(mesh.cornerIndices.size());
    d.cornerVertexIndices = mesh.cornerIndices.data();
    d.cornerWeights = mesh.cornerSharpnesses.data();
    d.numHoles = int(mesh.holeIndices.size());
    d.holeIndices = mesh.holeIndices.data();
    const bool hasTexcoords = mesh.texcoords.size() > 0 && mesh.texcoordIndices.size() == mesh.indices.size();
    Far::TopologyDescriptor::FVarChannel fvarChannel = {};
    if (hasTexcoords)
    {
        fvarChannel.numValues = int(mesh.texcoords.size());
        fvarChannel.valueIndices = mesh.texcoordIndices.data();
        d.numFVarChannels = 1;
        d.fvarChannels = &fvarChannel;
    }

    Sdc::Options sdcOptions;
    switch (mesh.interpolationRule)
    {
        case BOUNDARY_INTERPOLATION_NONE:
            sdcOptions.SetVtxBoundaryInterpolation(Sdc::Options::VTX_BOUNDARY_NONE);
            break;
        case BOUNDARY_INTERPOLATION_EDGE:
            sdcOptions.SetVtxBoundaryInterpolation(Sdc::Options::VTX_BOUNDARY_EDGE_ONLY);
            break;
        case BOUNDARY_INTERPOLATION_EDGE_AND_CORNER:
        default:
            sdcOptions.SetVtxBoundaryInterpolation(Sdc::Options::VTX_BOUNDARY_EDGE_AND_CORNER);
            break;
    }
    switch (mesh.fvarLinearInterpolation)
    {
        case FVAR_LINEAR_NONE:
            sdcOptions.SetFVarLinearInterpolation(Sdc::Options::FVAR_LINEAR_NONE);
            break;
        case FVAR_LINEAR_CORNERS_ONLY:
            sdcOptions.SetFVarLinearInterpolation(Sdc::Options::FVAR_LINEAR_CORNERS_ONLY);
            break;
        case FVAR_LINEAR_CORNERS_PLUS1:
            sdcOptions.SetFVarLinearInterpolation(Sdc::Options::FVAR_LINEAR_CORNERS_PLUS1);
            break;
        case FVAR_LINEAR_CORNERS_PLUS2:
            sdcOptions.SetFVarLinearInterpolation(Sdc::Options::FVAR_LINEAR_CORNERS_PLUS2);
            break;
        case FVAR_LINEAR_BOUNDARIES:
            sdcOptions.SetFVarLinearInterpolation(Sdc::Options::FVAR_LINEAR_BOUNDARIES);
            break;
        case FVAR_LINEAR_ALL:
            sdcOptions.SetFVarLinearInterpolation(Sdc::Options::FVAR_LINEAR_ALL);
            break;
        default:
            sdcOptions.SetFVarLinearInterpolation(Sdc::Options::FVAR_LINEAR_CORNERS_PLUS1);
            break;
    }
    sdcOptions.SetCreasingMethod(ToSdcCreasingMethod(options.creasingMethod));
    sdcOptions.SetTriangleSubdivision(ToSdcTriangleSubdivision(mesh.triangleSubdivisionRule));

    Far::TopologyRefinerFactory<Far::TopologyDescriptor>::Options o(
        SchemeFromString(options.subdivisionScheme), sdcOptions);

    SubdivisionEdgeMap edgeMap;
    edgeMap.reserve(mesh.indices.size());
    const int edgesWithOver2Faces = CountNonManifoldEdges(edgeMap);

    Far::TopologyRefiner *refiner =
        Far::TopologyRefinerFactory<Far::TopologyDescriptor>::Create(d, o);
    if (!refiner)
    {
        return false;
    }
    const bool hasRefinerFVar = refiner->GetNumFVarChannels() > 0;
    YBI_ERROR(!hasTexcoords || hasRefinerFVar,
              "Subdivision texcoords present but refiner has no fvar channels\n");

    int nextGeneratedVertexId = int(mesh.vertices.size());
    const std::vector<SubdivisionPatch> patches =
        BuildSubdivisionPatches(mesh, *refiner, edgeMap, nextGeneratedVertexId);

    const EdgeMapChecks edgeChecks = RunEdgeMapChecks(mesh, patches, edgeMap);
    if (!edgeChecks.ok)
    {
        delete refiner;
        return false;
    }

    const bool enableFVarTables = hasTexcoords && hasRefinerFVar;

    Far::TopologyRefiner::AdaptiveOptions adaptiveOptions(options.level);
    adaptiveOptions.considerFVarChannels = enableFVarTables;
    refiner->RefineAdaptive(adaptiveOptions);

    Far::PatchTableFactory::Options patchOptions(options.level);
    patchOptions.endCapType = Far::PatchTableFactory::Options::ENDCAP_GREGORY_BASIS;
    patchOptions.useInfSharpPatch = (d.numCreases > 0) || (d.numCorners > 0);
    patchOptions.generateFVarTables = enableFVarTables;
    patchOptions.generateFVarLegacyLinearPatches = false;
    const Far::PatchTable *patchTable = Far::PatchTableFactory::Create(*refiner, patchOptions);
    if (!patchTable)
    {
        delete refiner;
        return false;
    }

    Far::PatchMap patchMap(*patchTable);
    const std::vector<LimitEvalVertex> limitValues =
        BuildLimitEvalVertices(*refiner, *patchTable, mesh.vertices);
    std::vector<LimitEvalFVar2> limitFVarValues;
    const bool hasFVarChannel = hasTexcoords && hasRefinerFVar && patchTable->GetNumFVarChannels() > 0;
    if (hasFVarChannel)
    {
        limitFVarValues = BuildLimitEvalFVarValues(*refiner, *patchTable, mesh.texcoords, 0);
    }
    const bool enableLimitFVarEval = hasFVarChannel && !limitFVarValues.empty();

    const float3 eye = options.eye;
    const float3 lookAt = options.lookAt;

    int tmaxComputedEdges = 0;
    const std::vector<SubdivisionPatch> splitPatches =
        DiagSplitPatches(patches,
                         edgeMap,
                         nextGeneratedVertexId,
                         options.maxDiagSplitDepth,
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
                         options.useCameraMatrices,
                         options.cameraFromWorld,
                         options.clipFromCamera,
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
    int innerGridOnlyTris = 0;
    int stitchingOnlyTris = 0;
    std::vector<int> *triPatchFaceIds = nullptr;
    std::vector<int> *triCoarseFaceIds = nullptr;
    std::vector<int> *triPtexFaceIds = nullptr;
    std::vector<int> *triQuadrants = nullptr;
    if (options.generateTriangleMetadata)
    {
        triPatchFaceIds = &outResult->trianglePatchFaceIds;
        triCoarseFaceIds = &outResult->triangleCoarseFaceIds;
        triPtexFaceIds = &outResult->trianglePtexFaceIds;
        triQuadrants = &outResult->triangleQuadrants;
    }
    else
    {
        outResult->trianglePatchFaceIds.clear();
        outResult->triangleCoarseFaceIds.clear();
        outResult->trianglePtexFaceIds.clear();
        outResult->triangleQuadrants.clear();
    }
    if (!BuildLeafPatchStitchedMesh(splitPatches,
                                    edgeMap,
                                    patchMap,
                                    *patchTable,
                                    limitValues,
                                    enableLimitFVarEval ? &limitFVarValues : nullptr,
                                    nextGeneratedVertexId,
                                    &outResult->mesh,
                                    triPatchFaceIds,
                                    triCoarseFaceIds,
                                    triPtexFaceIds,
                                    triQuadrants,
                                    &innerGridVerts,
                                    &innerGridTris,
                                    &innerGridOnlyTris,
                                    &stitchingOnlyTris))
    {
        delete patchTable;
        delete refiner;
        return false;
    }

    outResult->refinedMaxLevel = refiner->GetMaxLevel();
    outResult->totalPatches = patchTable->GetNumPatchesTotal();
    outResult->subdivisionPatchCount = patches.size();
    outResult->diagSplitPatchCount = splitPatches.size();
    outResult->generatedVertexCount = nextGeneratedVertexId - int(mesh.vertices.size());
    outResult->midpointEdges = CountEdgesWithMidpointVertex(edgeMap);
    outResult->edgeTMaxComputed = tmaxComputedEdges;
    outResult->totalComputedEdges = CountEdgesWithComputedTMax(edgeMap);
    outResult->patchQuadVerts = patchQuadVerts;
    outResult->patchQuadCount = patchQuadCount;
    outResult->innerGridTriangleCount = innerGridOnlyTris;
    outResult->stitchingTriangleCount = stitchingOnlyTris;
    outResult->controlCageUniqueEdges = edgeMap.size();
    outResult->boundaryEdges = CountBoundaryEdges(edgeMap);
    outResult->controlCageEdgesWithOver2Faces = edgesWithOver2Faces;

    delete patchTable;
    delete refiner;
    return true;
}

YBI_NAMESPACE_END
