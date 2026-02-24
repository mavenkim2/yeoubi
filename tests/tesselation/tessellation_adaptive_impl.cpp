#include "io/usd_subdiv_json_io.h"
#include "tesselation/edge_map_validation.h"
#include "tesselation/subdivision_patch_types.h"

#include <opensubdiv/far/patchTableFactory.h>
#include <opensubdiv/far/ptexIndices.h>
#include <opensubdiv/far/topologyDescriptor.h>
#include <opensubdiv/far/topologyRefinerFactory.h>
#include <pxr/base/gf/vec2f.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <string>
#include <unordered_map>
#include <vector>

using namespace OpenSubdiv;

struct TriMesh
{
    pxr::VtVec3fArray positions;
    std::vector<int> indices;
};

static int ComputeDiagSplitPatchEdgeTMax(const std::vector<pxr::GfVec2f> &edgeSamples,
                                         float targetPixelSpacing)
{
    if (edgeSamples.size() < 2 || targetPixelSpacing <= 0.0f)
    {
        return 1;
    }

    float sumL = 0.0f;
    float maxL = 0.0f;
    for (size_t i = 1; i < edgeSamples.size(); ++i)
    {
        const pxr::GfVec2f d = edgeSamples[i] - edgeSamples[i - 1];
        const float l = std::sqrt(d[0] * d[0] + d[1] * d[1]);
        const float clamped = std::max(0.0f, l);
        sumL += clamped;
        maxL = std::max(maxL, clamped);
    }
    const int sampleCount = int(edgeSamples.size());
    const int tMax =
        std::max(1, int(std::ceil((float(sampleCount) * maxL) / targetPixelSpacing)));
    return tMax;
}

static int ComputePatchEdgeTMaxFactorsBasic(const std::vector<SubdivisionPatch> &patches,
                                            SubdivisionEdgeMap &edgeMap,
                                            float targetPixelSpacing)
{
    // Local patch UV domain starts at unit quad corners.
    const pxr::GfVec2f patchUVs[4] = {
        pxr::GfVec2f(0.0f, 0.0f),
        pxr::GfVec2f(1.0f, 0.0f),
        pxr::GfVec2f(1.0f, 1.0f),
        pxr::GfVec2f(0.0f, 1.0f)};

    int computedCount = 0;
    for (const SubdivisionPatch &patch : patches)
    {
        for (int edgeIndex = 0; edgeIndex < 4; ++edgeIndex)
        {
            const int next = (edgeIndex + 1) & 3;
            const uint64_t key = MakeEdgeKey(patch.verts[edgeIndex], patch.verts[next]);
            auto it = edgeMap.find(key);
            if (it == edgeMap.end())
            {
                continue;
            }

            SubdivisionEdge &edge = it->second;
            if (edge.tmaxComputed)
            {
                continue;
            }

            const std::vector<pxr::GfVec2f> samples = {patchUVs[edgeIndex], patchUVs[next]};
            edge.tmaxEdgeFactor = ComputeDiagSplitPatchEdgeTMax(samples, targetPixelSpacing);
            edge.tmaxComputed = true;
            computedCount++;
        }
    }
    return computedCount;
}

static SubdivisionEdgeMap BuildSubdivisionEdgeMap(const SelectedSubdivMesh &m)
{
    SubdivisionEdgeMap edgeMap;
    edgeMap.reserve(m.faceVertexIndices.size());

    size_t cursor = 0;
    for (size_t f = 0; f < m.faceVertexCounts.size(); ++f)
    {
        const int n = m.faceVertexCounts[f];
        if (n < 2 || cursor + n > m.faceVertexIndices.size())
        {
            cursor += std::max(0, n);
            continue;
        }

        for (int i = 0; i < n; ++i)
        {
            const int a = m.faceVertexIndices[cursor + i];
            const int b = m.faceVertexIndices[cursor + ((i + 1) % n)];
            const uint64_t key = MakeEdgeKey(a, b);
            SubdivisionEdge &edge = edgeMap[key];
            if (edge.faceCount == 0)
            {
                edge.v0 = std::min(a, b);
                edge.v1 = std::max(a, b);
                edge.firstFace = int(f);
            }
            else if (edge.faceCount == 1)
            {
                edge.secondFace = int(f);
            }
            edge.faceCount++;
        }
        cursor += n;
    }

    for (auto &kv : edgeMap)
    {
        SubdivisionEdge &edge = kv.second;
        edge.boundary = (edge.faceCount == 1);
        edge.nonManifold = (edge.faceCount > 2);
    }

    return edgeMap;
}

static SubdivisionEdge &GetOrCreateEdge(SubdivisionEdgeMap &edgeMap, int v0, int v1)
{
    const uint64_t key = MakeEdgeKey(v0, v1);
    SubdivisionEdge &edge = edgeMap[key];
    if (edge.v0 < 0 || edge.v1 < 0)
    {
        edge.v0 = std::min(v0, v1);
        edge.v1 = std::max(v0, v1);
    }
    return edge;
}

static int GetOrAllocateMidpointVertex(SubdivisionEdgeMap &edgeMap,
                                       int v0,
                                       int v1,
                                       int &nextVertexId)
{
    SubdivisionEdge &edge = GetOrCreateEdge(edgeMap, v0, v1);
    if (edge.midpointVertex < 0)
    {
        edge.midpointVertex = nextVertexId++;
    }
    return edge.midpointVertex;
}

static void AddEdgeUse(SubdivisionEdgeMap &edgeMap, int v0, int v1, int faceId)
{
    SubdivisionEdge &edge = GetOrCreateEdge(edgeMap, v0, v1);
    if (edge.faceCount == 0)
    {
        edge.firstFace = faceId;
    }
    else if (edge.faceCount == 1)
    {
        edge.secondFace = faceId;
    }
    edge.faceCount++;
}

static void FinalizeEdgeFlags(SubdivisionEdgeMap &edgeMap)
{
    for (auto &kv : edgeMap)
    {
        SubdivisionEdge &edge = kv.second;
        edge.boundary = (edge.faceCount == 1);
        edge.nonManifold = (edge.faceCount > 2);
    }
}

static std::vector<SubdivisionPatch> BuildSubdivisionPatches(const SelectedSubdivMesh &m,
                                                             const Far::TopologyRefiner &refiner,
                                                             SubdivisionEdgeMap &edgeMap,
                                                             int &nextVertexId)
{
    std::vector<SubdivisionPatch> patches;
    patches.reserve(m.faceVertexCounts.size());

    Far::PtexIndices ptexIndices(refiner);
    std::vector<int> faceCenterVertex(m.faceVertexCounts.size(), -1);

    size_t cursor = 0;
    for (size_t f = 0; f < m.faceVertexCounts.size(); ++f)
    {
        const int n = m.faceVertexCounts[f];
        if (n < 3 || cursor + n > m.faceVertexIndices.size())
        {
            cursor += std::max(0, n);
            continue;
        }

        const int basePtexFaceId = ptexIndices.GetFaceId(int(f));
        if (n == 4)
        {
            SubdivisionPatch patch = {};
            patch.verts[0] = m.faceVertexIndices[cursor + 0];
            patch.verts[1] = m.faceVertexIndices[cursor + 1];
            patch.verts[2] = m.faceVertexIndices[cursor + 2];
            patch.verts[3] = m.faceVertexIndices[cursor + 3];
            patch.coarseFace = int(f);
            patch.quadrant = 0;
            patch.ptexFaceId = basePtexFaceId;
            patches.push_back(patch);
        }
        else
        {
            if (faceCenterVertex[f] < 0)
            {
                faceCenterVertex[f] = nextVertexId++;
            }
            const int center = faceCenterVertex[f];

            for (int i = 0; i < n; ++i)
            {
                const int v = m.faceVertexIndices[cursor + i];
                const int vNext = m.faceVertexIndices[cursor + ((i + 1) % n)];
                const int vPrev = m.faceVertexIndices[cursor + ((i + n - 1) % n)];
                const int midNext = GetOrAllocateMidpointVertex(edgeMap, v, vNext, nextVertexId);
                const int midPrev = GetOrAllocateMidpointVertex(edgeMap, vPrev, v, nextVertexId);

                SubdivisionPatch patch = {};
                patch.verts[0] = v;
                patch.verts[1] = midNext;
                patch.verts[2] = center;
                patch.verts[3] = midPrev;
                patch.coarseFace = int(f);
                patch.quadrant = i;
                patch.ptexFaceId = basePtexFaceId + i;
                const int patchId = int(patches.size());
                patches.push_back(patch);

                // Only generated edges for non-quad quadrangulation.
                AddEdgeUse(edgeMap, v, midNext, patchId);
                AddEdgeUse(edgeMap, midNext, center, patchId);
                AddEdgeUse(edgeMap, center, midPrev, patchId);
                AddEdgeUse(edgeMap, midPrev, v, patchId);
            }
        }

        cursor += n;
    }

    FinalizeEdgeFlags(edgeMap);
    return patches;
}

struct CreasePairs
{
    std::vector<int> pairs;
    std::vector<float> weights;
};

static CreasePairs BuildCreasePairs(const SelectedSubdivMesh &m)
{
    CreasePairs out = {};
    size_t c = 0;
    for (size_t i = 0; i < m.creaseLengths.size(); i++)
    {
        const int len = m.creaseLengths[i];
        if (len < 2 || c + len > m.creaseIndices.size())
        {
            c += std::max(0, len);
            continue;
        }
        const float w = (i < m.creaseSharpnesses.size()) ? m.creaseSharpnesses[i] : 0.0f;
        for (int j = 0; j + 1 < len; j++)
        {
            out.pairs.push_back(m.creaseIndices[c + j]);
            out.pairs.push_back(m.creaseIndices[c + j + 1]);
            out.weights.push_back(w);
        }
        c += len;
    }
    return out;
}

static Sdc::SchemeType SchemeFromString(const std::string &s)
{
    if (s == "loop")
    {
        return Sdc::SCHEME_LOOP;
    }
    if (s == "bilinear")
    {
        return Sdc::SCHEME_BILINEAR;
    }
    return Sdc::SCHEME_CATMARK;
}

static Sdc::Options::VtxBoundaryInterpolation VtxBoundaryFromString(const std::string &s)
{
    if (s == "none")
    {
        return Sdc::Options::VTX_BOUNDARY_NONE;
    }
    if (s == "edgeOnly")
    {
        return Sdc::Options::VTX_BOUNDARY_EDGE_ONLY;
    }
    return Sdc::Options::VTX_BOUNDARY_EDGE_AND_CORNER;
}

static Sdc::Options::FVarLinearInterpolation FVarLinearFromString(const std::string &s)
{
    if (s == "none")
    {
        return Sdc::Options::FVAR_LINEAR_NONE;
    }
    if (s == "cornersOnly")
    {
        return Sdc::Options::FVAR_LINEAR_CORNERS_ONLY;
    }
    if (s == "cornersPlus1")
    {
        return Sdc::Options::FVAR_LINEAR_CORNERS_PLUS1;
    }
    if (s == "cornersPlus2")
    {
        return Sdc::Options::FVAR_LINEAR_CORNERS_PLUS2;
    }
    if (s == "boundaries")
    {
        return Sdc::Options::FVAR_LINEAR_BOUNDARIES;
    }
    if (s == "all" || s == "bilinear")
    {
        return Sdc::Options::FVAR_LINEAR_ALL;
    }
    return Sdc::Options::FVAR_LINEAR_CORNERS_PLUS1;
}

static Sdc::Options::CreasingMethod CreasingMethodFromString(const std::string &s)
{
    if (s == "chaikin")
    {
        return Sdc::Options::CREASE_CHAIKIN;
    }
    return Sdc::Options::CREASE_UNIFORM;
}

static Sdc::Options::TriangleSubdivision TriangleSubFromString(const std::string &s)
{
    if (s == "smooth")
    {
        return Sdc::Options::TRI_SUB_SMOOTH;
    }
    return Sdc::Options::TRI_SUB_CATMARK;
}

static bool WriteTriObj(const std::string &path, const TriMesh &mesh)
{
    std::ofstream out(path);
    if (!out.is_open())
    {
        return false;
    }
    for (const pxr::GfVec3f &p : mesh.positions)
    {
        out << "v " << p[0] << " " << p[1] << " " << p[2] << "\n";
    }
    for (int i = 0; i + 2 < int(mesh.indices.size()); i += 3)
    {
        out << "f " << (mesh.indices[i + 0] + 1) << " " << (mesh.indices[i + 1] + 1) << " "
            << (mesh.indices[i + 2] + 1) << "\n";
    }
    return true;
}

static TriMesh BuildTriangulatedControlCage(const SelectedSubdivMesh &m)
{
    TriMesh mesh = {};
    mesh.positions = m.points;

    size_t cursor = 0;
    for (size_t f = 0; f < m.faceVertexCounts.size(); ++f)
    {
        const int n = m.faceVertexCounts[f];
        if (n < 3 || cursor + n > m.faceVertexIndices.size())
        {
            cursor += std::max(0, n);
            continue;
        }

        const int i0 = m.faceVertexIndices[cursor + 0];
        for (int i = 1; i + 1 < n; ++i)
        {
            mesh.indices.push_back(i0);
            mesh.indices.push_back(m.faceVertexIndices[cursor + i]);
            mesh.indices.push_back(m.faceVertexIndices[cursor + i + 1]);
        }
        cursor += n;
    }

    return mesh;
}

static int CountNonManifoldEdges(const SubdivisionEdgeMap &edgeMap)
{
    int count = 0;
    for (const auto &it : edgeMap)
    {
        if (it.second.nonManifold)
        {
            count++;
        }
    }
    return count;
}

static int CountBoundaryEdges(const SubdivisionEdgeMap &edgeMap)
{
    int count = 0;
    for (const auto &it : edgeMap)
    {
        if (it.second.boundary)
        {
            count++;
        }
    }
    return count;
}

static int CountEdgesWithMidpointVertex(const SubdivisionEdgeMap &edgeMap)
{
    int count = 0;
    for (const auto &it : edgeMap)
    {
        if (it.second.midpointVertex >= 0)
        {
            count++;
        }
    }
    return count;
}

static int CountEdgesWithComputedTMax(const SubdivisionEdgeMap &edgeMap)
{
    int count = 0;
    for (const auto &it : edgeMap)
    {
        if (it.second.tmaxComputed)
        {
            count++;
        }
    }
    return count;
}

int main(int argc, char **argv)
{
    if (argc < 2)
    {
        std::fprintf(stderr, "Usage: %s <selected-subdiv.json> [level>=1] [out.obj]\n", argv[0]);
        return 2;
    }

    const std::string inJson = argv[1];
    int level = 1;
    std::string outObj;
    if (argc >= 3)
    {
        level = std::max(1, std::atoi(argv[2]));
    }
    if (argc >= 4)
    {
        outObj = argv[3];
    }
    if (outObj.empty())
    {
        outObj = "tests/bvh/out/refined_adaptive_level" + std::to_string(level) + ".obj";
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

    SubdivisionEdgeMap edgeMap = BuildSubdivisionEdgeMap(m);
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
    const int tmaxComputedEdges =
        ComputePatchEdgeTMaxFactorsBasic(patches, edgeMap, /*targetPixelSpacing*/ 0.25f);
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

    // TODO: restart adaptive tessellation path from patch params; current fallback writes
    // triangulated control cage while keeping refiner/patch table creation validated.
    const TriMesh triMesh = BuildTriangulatedControlCage(m);
    if (!WriteTriObj(outObj, triMesh))
    {
        std::fprintf(stderr, "Failed to write OBJ: %s\n", outObj.c_str());
        delete patchTable;
        delete refiner;
        return 1;
    }

    std::printf("Wrote fallback adaptive OBJ (control cage triangulation): %s\n", outObj.c_str());
    std::printf("  levelRequested=%d\n", level);
    std::printf("  fallback verts=%zu tris=%zu\n",
                triMesh.positions.size(),
                triMesh.indices.size() / 3);
    std::printf("  refinedMaxLevel=%d patches=%d\n",
                refiner->GetMaxLevel(),
                patchTable->GetNumPatchesTotal());
    std::printf("  subdivisionPatches=%zu generatedVertexCount=%d midpointEdges=%d\n",
                patches.size(),
                nextGeneratedVertexId - int(m.points.size()),
                CountEdgesWithMidpointVertex(edgeMap));
    std::printf("  edgeTMaxComputed=%d totalComputedEdges=%d\n",
                tmaxComputedEdges,
                CountEdgesWithComputedTMax(edgeMap));
    std::printf("  edgeMapChecks=ok\n");
    std::printf("  controlCageUniqueEdges=%zu boundaryEdges=%d\n",
                edgeMap.size(),
                CountBoundaryEdges(edgeMap));
    std::printf("  controlCageEdgesWithOver2Faces=%d\n", edgesWithOver2Faces);

    delete patchTable;
    delete refiner;
    return 0;
}
