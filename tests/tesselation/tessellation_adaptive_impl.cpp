#include "io/usd_subdiv_json_io.h"
#include "tesselation/edge_map_validation.h"
#include "tesselation/subdivision_patch_types.h"
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
#include <cstdlib>
#include <limits>
#include <string>
#include <unordered_map>
#include <vector>

using namespace OpenSubdiv;

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


#include "tesselation/tessellation_adaptive_limit_eval.inc"
#include "tesselation/tessellation_adaptive_edge_ops.inc"
#include "tesselation/tessellation_adaptive_patch_build.inc"
#include "tesselation/tessellation_adaptive_util.inc"
#include "tesselation/tessellation_adaptive_obj_write.inc"

int main(int argc, char **argv)
{
    constexpr float kDefaultPixelSpacing = 1.0f;
    constexpr int kDefaultSplitThreshold = 1;
    constexpr int kDefaultSampleSteps = 8;
    const std::string kDefaultPatchQuadObj = "tests/bvh/out/diagsplit_leaf_patch_quads.obj";
    const std::string kDefaultInnerGridObj = "tests/bvh/out/diagsplit_leaf_inner_grid.obj";

    if (argc < 2)
    {
        std::fprintf(stderr,
                     "Usage: %s <selected-subdiv.json> [level>=1]"
                     " [pixelSpacing>0] [splitThreshold>=1] [sampleSteps>=2] [outObjPath]\n",
                     argv[0]);
        return 2;
    }

    const std::string inJson = argv[1];
    int level = 1;
    float pixelSpacing = kDefaultPixelSpacing;
    int splitThreshold = kDefaultSplitThreshold;
    int sampleSteps = kDefaultSampleSteps;
    std::string outObjPath = kDefaultInnerGridObj;
    if (argc >= 3)
    {
        level = std::max(1, std::atoi(argv[2]));
    }
    if (argc >= 4)
    {
        pixelSpacing = std::max(1e-6f, float(std::atof(argv[3])));
    }
    if (argc >= 5)
    {
        splitThreshold = std::max(1, std::atoi(argv[4]));
    }
    if (argc >= 6)
    {
        sampleSteps = std::max(2, std::atoi(argv[5]));
    }
    if (argc >= 7)
    {
        outObjPath = argv[6];
    }

    SelectedSubdivMesh m = {};
    UsdCameraInfo camera = {};
    if (!ybi::testio::LoadSelectedSubdivFromJson(inJson, m, camera))
    {
        std::fprintf(stderr, "Failed to load JSON: %s\n", inJson.c_str());
        return 1;
    }

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

    // SubdivisionEdgeMap edgeMap = BuildSubdivisionEdgeMap(m);
    SubdivisionEdgeMap edgeMap;
    edgeMap.reserve(m.faceVertexIndices.size());
    const int edgesWithOver2Faces = CountNonManifoldEdges(edgeMap);
    if (edgesWithOver2Faces > 0)
    {
        std::fprintf(
            stderr,
            "WARNING: non-manifold control cage detected from indices (edges with >2 faces: %d)\n",
            edgesWithOver2Faces);
    }

    Far::TopologyRefiner *refiner =
        Far::TopologyRefinerFactory<Far::TopologyDescriptor>::Create(d, o);
    if (!refiner)
    {
        std::fprintf(stderr, "Failed to create TopologyRefiner\n");
        return 1;
    }

    int nextGeneratedVertexId = int(m.points.size());
    const std::vector<SubdivisionPatch> patches =
        BuildSubdivisionPatches(m, *refiner, edgeMap, nextGeneratedVertexId);
    const EdgeMapChecks edgeChecks = RunEdgeMapChecks(m, patches, edgeMap);
    if (!edgeChecks.ok)
    {
        std::fprintf(stderr,
                     "Edge map checks failed:"
                     " missingNgonMidpoints=%d"
                     " duplicateMidpointVertices=%d"
                     " missingGeneratedPatchEdges=%d"
                     " badBoundaryFlags=%d"
                     " badNonManifoldFlags=%d\n",
                     edgeChecks.missingNgonMidpoints,
                     edgeChecks.duplicateMidpointVertices,
                     edgeChecks.missingGeneratedPatchEdges,
                     edgeChecks.badBoundaryFlags,
                     edgeChecks.badNonManifoldFlags);
        delete refiner;
        return 1;
    }

    Far::TopologyRefiner::AdaptiveOptions adaptiveOptions(level);
    refiner->RefineAdaptive(adaptiveOptions);
    if (refiner->GetMaxLevel() < level)
    {
        std::fprintf(stderr, "Adaptive refine did not produce level %d\n", level);
        delete refiner;
        return 1;
    }

    Far::PatchTableFactory::Options patchOptions(level);
    patchOptions.endCapType = Far::PatchTableFactory::Options::ENDCAP_GREGORY_BASIS;
    patchOptions.useInfSharpPatch = (d.numCreases > 0) || (d.numCorners > 0);
    patchOptions.generateFVarLegacyLinearPatches = false;
    const Far::PatchTable *patchTable = Far::PatchTableFactory::Create(*refiner, patchOptions);
    if (!patchTable)
    {
        std::fprintf(stderr, "Failed to create PatchTable\n");
        delete refiner;
        return 1;
    }

    Far::PatchMap patchMap(*patchTable);
    const std::vector<LimitEvalVertex> limitValues =
        BuildLimitEvalVertices(*refiner, *patchTable, m.points);

    const pxr::GfVec3f meshCenter = ComputeMeshCenter(m.points);
    const pxr::GfVec3f eye =
        camera.found ? camera.worldPosition : (meshCenter + pxr::GfVec3f(0.0f, 0.0f, 5.0f));
    const pxr::GfVec3f lookAt = camera.found ? camera.meshCenter : meshCenter;
    // int tmaxComputedEdges = PrecomputePatchEdgeFactors(m,
    //                                                    patches,
    //                                                    edgeMap,
    //                                                    nextGeneratedVertexId,
    //                                                    patchMap,
    //                                                    *patchTable,
    //                                                    limitValues,
    //                                                    sampleSteps,
    //                                                    pixelSpacing,
    //                                                    splitThreshold,
    //                                                    eye,
    //                                                    lookAt,
    //                                                    /*viewportWidth*/ 1920,
    //                                                    /*viewportHeight*/ 1080,
    //                                                    /*verticalFovDegrees*/ 45.0f);
    int tmaxComputedEdges = 0;
    const std::vector<SubdivisionPatch> splitPatches =
        DiagSplitPatches(m,
                         patches,
                         edgeMap,
                         nextGeneratedVertexId,
                         patchMap,
                         *patchTable,
                         limitValues,
                         sampleSteps,
                         pixelSpacing,
                         splitThreshold,
                         eye,
                         lookAt,
                         /*viewportWidth*/ 1920,
                         /*viewportHeight*/ 1080,
                         /*verticalFovDegrees*/ 45.0f,
                         &tmaxComputedEdges);

    int missingStoredPatchParams = 0;
    int badUniformFactor = 0;
    if (!VerifyInitializedUniformEdgesHaveStoredPatchParams(
            edgeMap, &missingStoredPatchParams, &badUniformFactor))
    {
        std::fprintf(stderr,
                     "Initialized-uniform edge validation failed:"
                     " missingStoredPatchParams=%d"
                     " badUniformFactor=%d\n",
                     missingStoredPatchParams,
                     badUniformFactor);
        delete patchTable;
        delete refiner;
        return 1;
    }

    int innerGridVerts = 0;
    int innerGridTris = 0;
    int patchQuadVerts = 0;
    int patchQuadCount = 0;
    if (!WriteLeafPatchCornerQuadsObj(splitPatches,
                                      edgeMap,
                                      patchMap,
                                      *patchTable,
                                      limitValues,
                                      kDefaultPatchQuadObj,
                                      &patchQuadVerts,
                                      &patchQuadCount))
    {
        std::fprintf(
            stderr, "Failed to write leaf patch-quad OBJ: %s\n", kDefaultPatchQuadObj.c_str());
        delete patchTable;
        delete refiner;
        return 1;
    }
    if (!WriteLeafPatchInnerGridObj(splitPatches,
                                    edgeMap,
                                    patchMap,
                                    *patchTable,
                                    limitValues,
                                    outObjPath,
                                    &innerGridVerts,
                                    &innerGridTris))
    {
        std::fprintf(stderr, "Failed to write leaf inner-grid OBJ: %s\n", outObjPath.c_str());
        delete patchTable;
        delete refiner;
        return 1;
    }

    std::printf("  levelRequested=%d\n", level);
    std::printf("  refinedMaxLevel=%d patches=%d\n",
                refiner->GetMaxLevel(),
                patchTable->GetNumPatchesTotal());
    std::printf(
        "  subdivisionPatches=%zu diagSplitPatches=%zu generatedVertexCount=%d midpointEdges=%d\n",
        patches.size(),
        splitPatches.size(),
        nextGeneratedVertexId - int(m.points.size()),
        CountEdgesWithMidpointVertex(edgeMap));
    std::printf("  edgeTMaxComputed=%d totalComputedEdges=%d\n",
                tmaxComputedEdges,
                CountEdgesWithComputedTMax(edgeMap));
    std::printf("  diagSplitConfig pixelSpacing=%.3f splitThreshold=%d sampleSteps=%d\n",
                pixelSpacing,
                splitThreshold,
                sampleSteps);
    std::printf("  patchQuadObj path=%s verts=%d quads=%d\n",
                kDefaultPatchQuadObj.c_str(),
                patchQuadVerts,
                patchQuadCount);
    std::printf("  innerGridObj path=%s verts=%d tris=%d\n",
                outObjPath.c_str(),
                innerGridVerts,
                innerGridTris);
    std::printf("  edgeMapChecks=ok\n");
    std::printf("  controlCageUniqueEdges=%zu boundaryEdges=%d\n",
                edgeMap.size(),
                CountBoundaryEdges(edgeMap));
    std::printf("  controlCageEdgesWithOver2Faces=%d\n", edgesWithOver2Faces);

    delete patchTable;
    delete refiner;
    return 0;
}
