#include "io/usd_subdiv_json_io.h"

#include <opensubdiv/far/patchMap.h>
#include <opensubdiv/far/patchTable.h>
#include <opensubdiv/far/patchTableFactory.h>
#include <opensubdiv/far/primvarRefiner.h>
#include <opensubdiv/far/ptexIndices.h>
#include <opensubdiv/far/topologyDescriptor.h>
#include <opensubdiv/far/topologyRefinerFactory.h>

#include "util/assert.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <string>
#include <unordered_map>
#include <vector>

using namespace OpenSubdiv;

struct CreasePairs
{
    std::vector<int> pairs;
    std::vector<float> weights;
};

struct Primvar3
{
    float x = 0.0f;
    float y = 0.0f;
    float z = 0.0f;
    Primvar3() = default;
    explicit Primvar3(const pxr::GfVec3f &p) : x(p[0]), y(p[1]), z(p[2]) {}
    void Clear()
    {
        x = y = z = 0.0f;
    }
    void AddWithWeight(const Primvar3 &p, float w)
    {
        x += p.x * w;
        y += p.y * w;
        z += p.z * w;
    }
    pxr::GfVec3f ToGf() const
    {
        return pxr::GfVec3f(x, y, z);
    }
};

struct PtexFaceAdj
{
    std::array<int, 4> adjFace = {-1, -1, -1, -1};
    std::array<int, 4> adjEdge = {-1, -1, -1, -1};
    std::array<int, 4> edgeIndex = {-1, -1, -1, -1};
};

struct EdgeFactorCamera
{
    pxr::GfVec3f origin = pxr::GfVec3f(0.0f);
    pxr::GfVec3f target = pxr::GfVec3f(0.0f);
    pxr::GfVec3f up = pxr::GfVec3f(0.0f, 1.0f, 0.0f);
    float fovYDegrees = 50.0f;
    int width = 1920;
    int height = 1080;
};

struct EdgeFactorSettings
{
    int numSamples = 4;
    float pixelSpacing = 1.0f;
    int minRate = 1;
    int maxRate = 32;
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

static bool WriteObjLevel(const std::string &path,
                          const Far::TopologyLevel &level,
                          const std::vector<Primvar3> &pos,
                          int &nonQuadFacesOut)
{
    std::ofstream out(path);
    if (!out.is_open())
    {
        return false;
    }
    for (const auto &p : pos)
    {
        const pxr::GfVec3f g = p.ToGf();
        out << "v " << g[0] << " " << g[1] << " " << g[2] << "\n";
    }
    nonQuadFacesOut = 0;
    const int numVerts = int(pos.size());
    for (int f = 0; f < level.GetNumFaces(); f++)
    {
        auto verts = level.GetFaceVertices(f);
        if (verts.size() != 4)
        {
            nonQuadFacesOut++;
        }
        for (int i = 0; i < verts.size(); i++)
        {
            if (verts[i] < 0 || verts[i] >= numVerts)
            {
                std::fprintf(stderr,
                             "Invalid face index: face=%d corner=%d index=%d numVerts=%d\n",
                             f,
                             i,
                             verts[i],
                             numVerts);
                return false;
            }
        }
        out << "f";
        for (int i = 0; i < verts.size(); i++)
        {
            out << " " << (verts[i] + 1);
        }
        out << "\n";
    }
    return true;
}

static int CountControlCageEdgesWithOver2FacesFromIndices(const SelectedSubdivMesh &m)
{
    std::unordered_map<uint64_t, int> edgeUseCount;
    edgeUseCount.reserve(m.faceVertexIndices.size());

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
            const int lo = std::min(a, b);
            const int hi = std::max(a, b);
            const uint64_t key = (uint64_t(uint32_t(lo)) << 32) | uint64_t(uint32_t(hi));
            edgeUseCount[key]++;
        }
        cursor += n;
    }

    int over2 = 0;
    for (const auto &it : edgeUseCount)
    {
        if (it.second > 2)
        {
            over2++;
        }
    }
    return over2;
}

static pxr::GfVec3f Normalize(const pxr::GfVec3f &v)
{
    const float len = v.GetLength();
    if (len <= 1e-8f)
    {
        return pxr::GfVec3f(0.0f, 0.0f, 1.0f);
    }
    return v / len;
}

static EdgeFactorCamera BuildEdgeFactorCamera(const pxr::VtVec3fArray &points,
                                              const UsdCameraInfo &cameraInfo)
{
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

    EdgeFactorCamera camera = {};
    camera.target = (bbMin + bbMax) * 0.5f;
    const float diag = std::max((bbMax - bbMin).GetLength(), 1.0f);
    const float distance = cameraInfo.found ? std::max(cameraInfo.distanceToMeshCenter, 1e-3f)
                                            : (2.8f * diag);
    const pxr::GfVec3f dir = Normalize(pxr::GfVec3f(0.0f, 0.2f, 2.8f));
    camera.origin = camera.target + dir * distance;
    return camera;
}

static pxr::GfVec2f ToScreen(const EdgeFactorCamera &camera, const pxr::GfVec3f &point)
{
    const pxr::GfVec3f forward = Normalize(camera.target - camera.origin);
    const pxr::GfVec3f right = Normalize(pxr::GfCross(forward, camera.up));
    const pxr::GfVec3f up = Normalize(pxr::GfCross(right, forward));

    const pxr::GfVec3f rel = point - camera.origin;
    const float x = pxr::GfDot(rel, right);
    const float y = pxr::GfDot(rel, up);
    const float z = std::max(pxr::GfDot(rel, forward), 1e-4f);

    const float aspect = float(camera.width) / float(camera.height);
    const float tanHalfFovY = tanf(0.5f * camera.fovYDegrees * 0.0174532925199432957692f);
    const float tanHalfFovX = tanHalfFovY * aspect;

    const float ndcX = x / (z * tanHalfFovX);
    const float ndcY = y / (z * tanHalfFovY);
    const float sx = (ndcX * 0.5f + 0.5f) * float(camera.width);
    const float sy = (1.0f - (ndcY * 0.5f + 0.5f)) * float(camera.height);
    return pxr::GfVec2f(sx, sy);
}

static std::vector<Primvar3> BuildPatchEvalPositions(const Far::TopologyRefiner &refiner,
                                                      const Far::PatchTable &patchTable,
                                                      const pxr::VtVec3fArray &points)
{
    const int numVertices = refiner.GetNumVerticesTotal();
    const int numLocalPoints = patchTable.GetNumLocalPoints();
    std::vector<Primvar3> values(size_t(numVertices + numLocalPoints));
    for (size_t i = 0; i < points.size(); ++i)
    {
        values[i] = Primvar3(points[i]);
    }
    Primvar3 *src = values.data();
    Far::PrimvarRefiner prim(refiner);
    for (int l = 1; l < refiner.GetNumLevels(); ++l)
    {
        Primvar3 *dst = src + refiner.GetLevel(l - 1).GetNumVertices();
        prim.Interpolate(l, src, dst);
        src = dst;
    }
    patchTable.ComputeLocalPointValues(values.data(), values.data() + numVertices);
    return values;
}

static pxr::GfVec2f EdgeUV(int edge, float t)
{
    switch (edge)
    {
        case 0:
            return pxr::GfVec2f(t, 0.0f);
        case 1:
            return pxr::GfVec2f(1.0f, t);
        case 2:
            return pxr::GfVec2f(1.0f - t, 1.0f);
        default:
            return pxr::GfVec2f(0.0f, 1.0f - t);
    }
}

static bool EvaluateLimitPosition(const Far::PatchMap &patchMap,
                                  const Far::PatchTable &patchTable,
                                  const std::vector<Primvar3> &positions,
                                  int ptexFace,
                                  float u,
                                  float v,
                                  pxr::GfVec3f &outPos)
{
    const Far::PatchTable::PatchHandle *handle = patchMap.FindPatch(ptexFace, u, v);
    if (!handle)
    {
        const float eps = 1e-5f;
        u = std::clamp(u, eps, 1.0f - eps);
        v = std::clamp(v, eps, 1.0f - eps);
        handle = patchMap.FindPatch(ptexFace, u, v);
        if (!handle)
        {
            return false;
        }
    }

    const Far::ConstIndexArray cvIndices = patchTable.GetPatchVertices(*handle);
    std::vector<float> pWeights(size_t(cvIndices.size()));
    patchTable.EvaluateBasis(*handle, u, v, pWeights.data());

    Primvar3 p = {};
    for (int cv = 0; cv < cvIndices.size(); ++cv)
    {
        p.AddWithWeight(positions[size_t(cvIndices[cv])], pWeights[size_t(cv)]);
    }
    YBI_ASSERT(p.x != 0.0f || p.y != 0.0f || p.z != 0.0f);
    outPos = p.ToGf();
    return true;
}

static int ComputeDiagSplitEdgeFactor(const Far::PatchMap &patchMap,
                                      const Far::PatchTable &patchTable,
                                      const std::vector<Primvar3> &positions,
                                      int ptexFace,
                                      int edge,
                                      const EdgeFactorCamera &camera,
                                      const EdgeFactorSettings &settings)
{
    const int sampleCount = std::max(2, settings.numSamples);
    float maxLi = 0.0f;
    float sumLi = 0.0f;

    const pxr::GfVec2f uv0 = EdgeUV(edge, 0.0f);
    pxr::GfVec3f p0(0.0f);
    if (!EvaluateLimitPosition(patchMap, patchTable, positions, ptexFace, uv0[0], uv0[1], p0))
    {
        return settings.minRate;
    }
    pxr::GfVec2f prev = ToScreen(camera, p0);
    for (int i = 1; i < sampleCount; ++i)
    {
        const float t = float(i) / float(sampleCount - 1);
        const pxr::GfVec2f uv = EdgeUV(edge, t);
        pxr::GfVec3f p(0.0f);
        if (!EvaluateLimitPosition(patchMap, patchTable, positions, ptexFace, uv[0], uv[1], p))
        {
            return settings.minRate;
        }
        const pxr::GfVec2f s = ToScreen(camera, p);
        const pxr::GfVec2f d = s - prev;
        const float li = sqrtf(d[0] * d[0] + d[1] * d[1]);
        sumLi += li;
        maxLi = std::max(maxLi, li);
        prev = s;
    }

    const int tMin = int(ceilf(sumLi / settings.pixelSpacing));
    const int tMax = int(ceilf(float(sampleCount - 1) * maxLi / settings.pixelSpacing));
    (void)tMin;
    return std::clamp(std::max(1, tMax), settings.minRate, settings.maxRate);
}

static int ComputeMaxPtexEdgeFactor(const Far::PatchMap &patchMap,
                                    const Far::PatchTable &patchTable,
                                    const std::vector<Primvar3> &positions,
                                    const std::vector<PtexFaceAdj> &faces,
                                    int numPtexFaces,
                                    int uniqueEdgeCount,
                                    const EdgeFactorCamera &camera,
                                    const EdgeFactorSettings &settings)
{
    std::vector<int> edgeFactors(std::max(0, uniqueEdgeCount), -1);

    int maxCalculatedEdgeFactor = 0;
    for (int pf = 0; pf < numPtexFaces; ++pf)
    {
        if (!patchMap.FindPatch(pf, 0.5f, 0.5f))
        {
            continue;
        }
        for (int e = 0; e < 4; ++e)
        {
            const int edgeId = faces[pf].edgeIndex[e];
            YBI_ASSERT(edgeId >= 0 && edgeId < uniqueEdgeCount);
            const int candidate =
                ComputeDiagSplitEdgeFactor(patchMap, patchTable, positions, pf, e, camera, settings);
            if (edgeFactors[edgeId] < 0)
            {
                edgeFactors[edgeId] = candidate;
            }
            else
            {
                edgeFactors[edgeId] = std::max(edgeFactors[edgeId], candidate);
            }
            maxCalculatedEdgeFactor = std::max(maxCalculatedEdgeFactor, edgeFactors[edgeId]);
        }
    }
    return maxCalculatedEdgeFactor;
}

static int BuildUniquePtexEdgeIds(const Far::TopologyRefiner &refiner,
                                  const Far::TopologyLevel &level0,
                                  const Far::PtexIndices &ptex,
                                  std::vector<PtexFaceAdj> &facesOut,
                                  int &numPtexFacesOut)
{
    numPtexFacesOut = ptex.GetNumFaces();
    facesOut.assign(std::max(0, numPtexFacesOut), PtexFaceAdj{});

    for (int f = 0; f < level0.GetNumFaces(); ++f)
    {
        const int n = level0.GetFaceVertices(f).size();
        const int qCount = (n == 4) ? 1 : n;
        const int basePtex = ptex.GetFaceId(f);
        for (int q = 0; q < qCount; ++q)
        {
            const int pf = basePtex + q;
            YBI_ASSERT(pf >= 0 && pf < numPtexFacesOut);
            int adjFaces[4] = {-1, -1, -1, -1};
            int adjEdges[4] = {-1, -1, -1, -1};
            ptex.GetAdjacency(refiner, f, q, adjFaces, adjEdges);
            for (int e = 0; e < 4; ++e)
            {
                facesOut[pf].adjFace[e] = adjFaces[e];
                facesOut[pf].adjEdge[e] = adjEdges[e];
            }
        }
    }

    int edgeCount = 0;
    for (int pf = 0; pf < numPtexFacesOut; ++pf)
    {
        for (int e = 0; e < 4; ++e)
        {
            const int af = facesOut[pf].adjFace[e];
            const int ae = facesOut[pf].adjEdge[e];

            if (facesOut[pf].edgeIndex[e] >= 0)
            {
                continue;
            }

            if (af == -1)
            {
                facesOut[pf].edgeIndex[e] = edgeCount++;
                continue;
            }

            YBI_ASSERT(af >= 0 && af < numPtexFacesOut && ae >= 0 && ae < 4);

            const int neighborEdge = facesOut[af].edgeIndex[ae];
            if (neighborEdge >= 0)
            {
                facesOut[pf].edgeIndex[e] = neighborEdge;
                continue;
            }
            const int edgeId = edgeCount++;
            facesOut[pf].edgeIndex[e] = edgeId;
            facesOut[af].edgeIndex[ae] = edgeId;
        }
    }
    return edgeCount;
}

int main(int argc, char **argv)
{
    if (argc < 2 || argc > 4)
    {
        std::fprintf(stderr, "Usage: %s <selected-subdiv.json> [level>=1] [out.obj]\n", argv[0]);
        return 2;
    }

    const std::string inJson = argv[1];
    const int level = (argc >= 3) ? std::max(1, std::atoi(argv[2])) : 1;
    const std::string outObj =
        (argc >= 4) ? argv[3]
                    : ("tests/bvh/out/refined_adaptive_level" + std::to_string(level) + ".obj");

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

    Far::TopologyRefinerFactory<Far::TopologyDescriptor>::Options o(
        SchemeFromString(m.subdivisionScheme));

    const int edgesWithOver2Faces = CountControlCageEdgesWithOver2FacesFromIndices(m);
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

    const Far::TopologyLevel &level0 = refiner->GetLevel(0);

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
    const std::vector<Primvar3> patchEvalPositions = BuildPatchEvalPositions(*refiner, *patchTable, m.points);
    const EdgeFactorCamera edgeFactorCamera = BuildEdgeFactorCamera(m.points, camera);
    const EdgeFactorSettings edgeFactorSettings = {};

    std::vector<Primvar3> p0;
    p0.reserve(m.points.size());
    for (const auto &p : m.points)
    {
        p0.emplace_back(p);
    }
    std::vector<Primvar3> p1(refiner->GetLevel(level).GetNumVertices());
    Far::PrimvarRefiner prim(*refiner);
    prim.Interpolate(level, p0, p1);

    Far::PtexIndices ptex(*refiner);
    int ptexFaceCount = 0;
    std::vector<PtexFaceAdj> ptexFaceAdj;
    const int uniquePtexEdges =
        BuildUniquePtexEdgeIds(*refiner, level0, ptex, ptexFaceAdj, ptexFaceCount);
    const int maxCalculatedEdgeFactor = ComputeMaxPtexEdgeFactor(patchMap,
                                                                 *patchTable,
                                                                 patchEvalPositions,
                                                                 ptexFaceAdj,
                                                                 ptexFaceCount,
                                                                 uniquePtexEdges,
                                                                 edgeFactorCamera,
                                                                 edgeFactorSettings);

    int nonQuadFaces = 0;
    if (!WriteObjLevel(outObj, refiner->GetLevel(level), p1, nonQuadFaces))
    {
        std::fprintf(stderr, "Failed to write OBJ: %s\n", outObj.c_str());
        delete patchTable;
        delete refiner;
        return 1;
    }

    std::printf("Wrote adaptive level-%d OBJ: %s\n", level, outObj.c_str());
    std::printf("  verts=%zu faces=%d nonQuads=%d\n",
                p1.size(),
                refiner->GetLevel(level).GetNumFaces(),
                nonQuadFaces);
    std::printf("  maxCalculatedEdgeFactor=%d\n", maxCalculatedEdgeFactor);
    std::printf("  ptexFaces=%d\n", ptexFaceCount);
    std::printf("  controlCageEdgesWithOver2Faces=%d\n", edgesWithOver2Faces);

    delete patchTable;
    delete refiner;
    return 0;
}
