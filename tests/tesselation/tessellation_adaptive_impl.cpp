#include "io/usd_subdiv_json_io.h"
#include "tesselation/edge_map_validation.h"
#include "tesselation/subdivision_patch_types.h"

#include <opensubdiv/far/patchTableFactory.h>
#include <opensubdiv/far/patchMap.h>
#include <opensubdiv/far/ptexIndices.h>
#include <opensubdiv/far/primvarRefiner.h>
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

static float Dot(const pxr::GfVec3f &a, const pxr::GfVec3f &b)
{
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

static pxr::GfVec3f Cross(const pxr::GfVec3f &a, const pxr::GfVec3f &b)
{
    return pxr::GfVec3f(
        a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]);
}

static pxr::GfVec3f Normalize(const pxr::GfVec3f &v)
{
    const float lenSq = Dot(v, v);
    if (lenSq <= 1e-12f)
    {
        return pxr::GfVec3f(0.0f, 0.0f, 1.0f);
    }
    return v * (1.0f / std::sqrt(lenSq));
}

static std::vector<LimitEvalVertex> BuildLimitEvalVertices(const Far::TopologyRefiner &refiner,
                                                           const Far::PatchTable &patchTable,
                                                           const pxr::VtVec3fArray &coarsePoints)
{
    const int numRefinerVerts = refiner.GetNumVerticesTotal();
    const int numLocalPoints = patchTable.GetNumLocalPoints();
    std::vector<LimitEvalVertex> values(numRefinerVerts + numLocalPoints);

    const int numCoarseVerts = refiner.GetLevel(0).GetNumVertices();
    const int copyCount = std::min(numCoarseVerts, int(coarsePoints.size()));
    for (int i = 0; i < copyCount; ++i)
    {
        values[i].p = coarsePoints[i];
    }

    Far::PrimvarRefiner primvarRefiner(refiner);
    LimitEvalVertex *src = values.data();
    for (int level = 1; level < refiner.GetNumLevels(); ++level)
    {
        LimitEvalVertex *dst = src + refiner.GetLevel(level - 1).GetNumVertices();
        primvarRefiner.Interpolate(level, src, dst);
        src = dst;
    }

    if (numLocalPoints > 0)
    {
        patchTable.ComputeLocalPointValues(values.data(), values.data() + numRefinerVerts);
    }
    return values;
}

static pxr::GfVec3f EvaluateLimitPosition(const Far::PatchMap &patchMap,
                                          const Far::PatchTable &patchTable,
                                          const std::vector<LimitEvalVertex> &limitValues,
                                          int ptexFaceId,
                                          const pxr::GfVec2f &uv)
{
    const Far::PatchTable::PatchHandle *handle = patchMap.FindPatch(ptexFaceId, uv[0], uv[1]);
    if (!handle)
    {
        return pxr::GfVec3f(0.0f);
    }

    float pWeights[20] = {0.0f};
    patchTable.EvaluateBasis(*handle, uv[0], uv[1], pWeights);
    Far::ConstIndexArray cvs = patchTable.GetPatchVertices(*handle);

    pxr::GfVec3f p(0.0f);
    for (int i = 0; i < cvs.size(); ++i)
    {
        p += limitValues[cvs[i]].p * pWeights[i];
    }
    return p;
}

static pxr::GfVec2f ProjectToScreen(const pxr::GfVec3f &p,
                                    const pxr::GfVec3f &eye,
                                    const pxr::GfVec3f &lookAt,
                                    int viewportWidth,
                                    int viewportHeight,
                                    float verticalFovDegrees)
{
    const pxr::GfVec3f forward = Normalize(lookAt - eye);
    pxr::GfVec3f worldUp(0.0f, 0.0f, 1.0f);
    if (std::abs(Dot(forward, worldUp)) > 0.999f)
    {
        worldUp = pxr::GfVec3f(0.0f, 1.0f, 0.0f);
    }
    const pxr::GfVec3f right = Normalize(Cross(forward, worldUp));
    const pxr::GfVec3f up = Normalize(Cross(right, forward));

    const pxr::GfVec3f v = p - eye;
    const float x = Dot(v, right);
    const float y = Dot(v, up);
    const float z = std::max(1e-6f, Dot(v, forward));

    const float fovY = verticalFovDegrees * 3.14159265358979323846f / 180.0f;
    const float tanHalfFovY = std::tan(0.5f * fovY);
    const float aspect = float(viewportWidth) / float(viewportHeight);
    const float ndcX = x / (z * tanHalfFovY * aspect);
    const float ndcY = y / (z * tanHalfFovY);

    const float sx = (ndcX * 0.5f + 0.5f) * float(viewportWidth);
    const float sy = (1.0f - (ndcY * 0.5f + 0.5f)) * float(viewportHeight);
    return pxr::GfVec2f(sx, sy);
}

static int ComputeDiagSplitPatchEdgeTMax(const Far::PatchMap &patchMap,
                                         const Far::PatchTable &patchTable,
                                         const std::vector<LimitEvalVertex> &limitValues,
                                         int ptexFaceId,
                                         const pxr::GfVec2f &uvStart,
                                         const pxr::GfVec2f &uvEnd,
                                         int sampleStepsN,
                                         float targetPixelSpacing,
                                         const pxr::GfVec3f &eye,
                                         const pxr::GfVec3f &lookAt,
                                         int viewportWidth,
                                         int viewportHeight,
                                         float verticalFovDegrees)
{
    if (sampleStepsN < 2 || targetPixelSpacing <= 0.0f)
    {
        return 1;
    }

    float maxLi = 0.0f;
    const pxr::GfVec3f p0 = EvaluateLimitPosition(patchMap, patchTable, limitValues, ptexFaceId, uvStart);
    pxr::GfVec2f prev = ProjectToScreen(
        p0, eye, lookAt, viewportWidth, viewportHeight, verticalFovDegrees);
    for (int i = 1; i < sampleStepsN; ++i)
    {
        const float t = float(i) / float(sampleStepsN - 1);
        const pxr::GfVec2f uv = uvStart * (1.0f - t) + uvEnd * t;
        const pxr::GfVec3f p = EvaluateLimitPosition(patchMap, patchTable, limitValues, ptexFaceId, uv);
        const pxr::GfVec2f s =
            ProjectToScreen(p, eye, lookAt, viewportWidth, viewportHeight, verticalFovDegrees);
        const pxr::GfVec2f d = s - prev;
        const float li = std::sqrt(d[0] * d[0] + d[1] * d[1]);
        maxLi = std::max(maxLi, li);
        prev = s;
    }

    // DiagSplit Figure 8: tmax = ceil((N-1) * max(Li) / R).
    const int tMax =
        std::max(1, int(std::ceil(float(sampleStepsN - 1) * maxLi / targetPixelSpacing)));
    return tMax;
}

static int ComputePatchEdgeTMaxFactors(const std::vector<SubdivisionPatch> &patches,
                                       SubdivisionEdgeMap &edgeMap,
                                       const Far::PatchMap &patchMap,
                                       const Far::PatchTable &patchTable,
                                       const std::vector<LimitEvalVertex> &limitValues,
                                       int sampleStepsN,
                                       float targetPixelSpacing,
                                       const pxr::GfVec3f &eye,
                                       const pxr::GfVec3f &lookAt,
                                       int viewportWidth,
                                       int viewportHeight,
                                       float verticalFovDegrees)
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

            edge.tmaxEdgeFactor = ComputeDiagSplitPatchEdgeTMax(patchMap,
                                                                patchTable,
                                                                limitValues,
                                                                patch.ptexFaceId,
                                                                patchUVs[edgeIndex],
                                                                patchUVs[next],
                                                                sampleStepsN,
                                                                targetPixelSpacing,
                                                                eye,
                                                                lookAt,
                                                                viewportWidth,
                                                                viewportHeight,
                                                                verticalFovDegrees);
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
    const pxr::GfVec3f eye = camera.found ? camera.worldPosition : (meshCenter + pxr::GfVec3f(0.0f, 0.0f, 5.0f));
    const pxr::GfVec3f lookAt = camera.found ? camera.meshCenter : meshCenter;
    const int tmaxComputedEdges = ComputePatchEdgeTMaxFactors(patches,
                                                              edgeMap,
                                                              patchMap,
                                                              *patchTable,
                                                              limitValues,
                                                              /*sampleStepsN*/ 8,
                                                              /*targetPixelSpacing*/ 0.25f,
                                                              eye,
                                                              lookAt,
                                                              /*viewportWidth*/ 1920,
                                                              /*viewportHeight*/ 1080,
                                                              /*verticalFovDegrees*/ 45.0f);

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
