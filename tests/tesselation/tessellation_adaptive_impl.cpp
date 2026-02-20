#include "io/usd_subdiv_json_io.h"

#include <opensubdiv/far/patchMap.h>
#include <opensubdiv/far/patchTable.h>
#include <opensubdiv/far/patchTableFactory.h>
#include <opensubdiv/far/primvarRefiner.h>
#include <opensubdiv/far/ptexIndices.h>
#include <opensubdiv/far/topologyDescriptor.h>
#include <opensubdiv/far/topologyRefinerFactory.h>

#include "tesselation/edge_rate_obj_io.h"
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
    bool fromNgon = false;
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

struct EdgeFactorResult
{
    int maxCalculatedEdgeFactor = 0;
    std::vector<int> edgeFactors;
};

struct RefinedPositions
{
    std::vector<Primvar3> values;
    std::vector<int> levelStarts;
    int numRefinerVerts = 0;
};

struct TriMesh
{
    std::vector<pxr::GfVec3f> positions;
    std::vector<int> indices;
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
                                              const UsdCameraInfo &cameraInfo,
                                              float cameraDistanceScale)
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
    const float distance =
        (cameraInfo.found ? std::max(cameraInfo.distanceToMeshCenter, 1e-3f) : (2.8f * diag)) *
        std::max(cameraDistanceScale, 1e-4f);
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

static RefinedPositions BuildRefinedPositions(const Far::TopologyRefiner &refiner,
                                              const Far::PatchTable &patchTable,
                                              const pxr::VtVec3fArray &points)
{
    RefinedPositions out = {};
    const int numLevels = refiner.GetNumLevels();
    out.levelStarts.assign(std::max(0, numLevels), 0);
    for (int l = 1; l < numLevels; ++l)
    {
        out.levelStarts[l] = out.levelStarts[l - 1] + refiner.GetLevel(l - 1).GetNumVertices();
    }

    out.numRefinerVerts = refiner.GetNumVerticesTotal();
    const int numLocalPoints = patchTable.GetNumLocalPoints();
    out.values.resize(size_t(out.numRefinerVerts + numLocalPoints));
    for (size_t i = 0; i < points.size(); ++i)
    {
        out.values[i] = Primvar3(points[i]);
    }

    Primvar3 *src = out.values.data();
    Far::PrimvarRefiner prim(refiner);
    for (int l = 0; l < refiner.GetMaxLevel(); ++l)
    {
        Primvar3 *dst = src + refiner.GetLevel(l).GetNumVertices();
        prim.Interpolate(l + 1, src, dst);
        src = dst;
    }
    if (numLocalPoints > 0)
    {
        patchTable.ComputeLocalPointValues(out.values.data(), out.values.data() + out.numRefinerVerts);
    }
    return out;
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

static EdgeFactorResult ComputeEdgeFactors(const Far::PatchMap &patchMap,
                                           const Far::PatchTable &patchTable,
                                           const std::vector<Primvar3> &positions,
                                           const std::vector<PtexFaceAdj> &faces,
                                           int numPtexFaces,
                                           int uniqueEdgeCount,
                                           const EdgeFactorCamera &camera,
                                           const EdgeFactorSettings &settings)
{
    EdgeFactorResult result = {};
    result.edgeFactors.assign(std::max(0, uniqueEdgeCount), -1);
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
            if (result.edgeFactors[edgeId] < 0)
            {
                result.edgeFactors[edgeId] = candidate;
            }
            else
            {
                result.edgeFactors[edgeId] = std::max(result.edgeFactors[edgeId], candidate);
            }
            result.maxCalculatedEdgeFactor =
                std::max(result.maxCalculatedEdgeFactor, result.edgeFactors[edgeId]);
        }
    }
    return result;
}

static std::vector<ybi::testio::EdgeRateDebugLine>
BuildEdgeRateDebugLines(const Far::PatchMap &patchMap,
                        const Far::PatchTable &patchTable,
                        const std::vector<Primvar3> &positions,
                        const std::vector<PtexFaceAdj> &faces,
                        int numPtexFaces,
                        const std::vector<int> &edgeFactors)
{
    std::vector<ybi::testio::EdgeRateDebugLine> lines;
    std::vector<bool> emitted(edgeFactors.size(), false);

    for (int pf = 0; pf < numPtexFaces; ++pf)
    {
        if (!patchMap.FindPatch(pf, 0.5f, 0.5f))
        {
            continue;
        }
        for (int e = 0; e < 4; ++e)
        {
            const int edgeId = faces[pf].edgeIndex[e];
            if (edgeId < 0 || edgeId >= int(edgeFactors.size()) || emitted[edgeId])
            {
                continue;
            }
            emitted[edgeId] = true;

            const int rate = std::max(1, edgeFactors[edgeId]);
            const int samples = std::max(2, rate + 1);
            ybi::testio::EdgeRateDebugLine line = {};
            line.rate = rate;
            line.points.reserve(samples);

            bool valid = true;
            for (int i = 0; i < samples; ++i)
            {
                const float t = float(i) / float(samples - 1);
                const pxr::GfVec2f uv = EdgeUV(e, t);
                pxr::GfVec3f p(0.0f);
                if (!EvaluateLimitPosition(
                        patchMap, patchTable, positions, pf, uv[0], uv[1], p))
                {
                    valid = false;
                    break;
                }
                line.points.push_back(p);
            }
            if (valid && line.points.size() > 1)
            {
                lines.push_back(std::move(line));
            }
        }
    }
    return lines;
}

static TriMesh TessellateFixedRateNoStitch(const Far::PatchMap &patchMap,
                                           const Far::PatchTable &patchTable,
                                           const std::vector<Primvar3> &positions,
                                           int numPtexFaces,
                                           int fixedEdgeRate,
                                           int &skippedPtexFacesOut)
{
    TriMesh out = {};
    skippedPtexFacesOut = 0;
    const int steps = std::max(1, fixedEdgeRate);
    const int side = steps + 1;

    for (int pf = 0; pf < numPtexFaces; ++pf)
    {
        if (!patchMap.FindPatch(pf, 0.5f, 0.5f))
        {
            skippedPtexFacesOut++;
            continue;
        }

        const int base = int(out.positions.size());
        bool ok = true;
        for (int v = 0; v <= steps && ok; ++v)
        {
            for (int u = 0; u <= steps; ++u)
            {
                const float fu = float(u) / float(steps);
                const float fv = float(v) / float(steps);
                pxr::GfVec3f p(0.0f);
                if (!EvaluateLimitPosition(patchMap, patchTable, positions, pf, fu, fv, p))
                {
                    ok = false;
                    break;
                }
                out.positions.push_back(p);
            }
        }
        if (!ok)
        {
            out.positions.resize(size_t(base));
            skippedPtexFacesOut++;
            continue;
        }

        auto Idx = [side, base](int u, int v) { return base + v * side + u; };
        for (int v = 0; v < steps; ++v)
        {
            for (int u = 0; u < steps; ++u)
            {
                const int i0 = Idx(u, v);
                const int i1 = Idx(u + 1, v);
                const int i2 = Idx(u + 1, v + 1);
                const int i3 = Idx(u, v + 1);
                out.indices.push_back(i0);
                out.indices.push_back(i1);
                out.indices.push_back(i2);
                out.indices.push_back(i0);
                out.indices.push_back(i2);
                out.indices.push_back(i3);
            }
        }
    }
    return out;
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
            facesOut[pf].fromNgon = (n != 4);
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

static bool WriteMissingPatchNeighbors(const Far::PatchMap &patchMap,
                                       const std::vector<PtexFaceAdj> &faces,
                                       int numPtexFaces,
                                       const std::string &outPath,
                                       int &missingCountOut)
{
    std::ofstream out(outPath);
    if (!out.is_open())
    {
        return false;
    }

    missingCountOut = 0;
    for (int pf = 0; pf < numPtexFaces; ++pf)
    {
        if (patchMap.FindPatch(pf, 0.5f, 0.5f))
        {
            continue;
        }
        out << "ptexFace=" << pf;
        for (int e = 0; e < 4; ++e)
        {
            out << " e" << e << ":(adjFace=" << faces[pf].adjFace[e]
                << ",adjEdge=" << faces[pf].adjEdge[e] << ",edgeId=" << faces[pf].edgeIndex[e]
                << ")";
        }
        out << " fromNgon=" << (faces[pf].fromNgon ? 1 : 0) << "\n";
        missingCountOut++;
    }
    return out.good();
}

int main(int argc, char **argv)
{
    if (argc < 2)
    {
        std::fprintf(stderr,
                     "Usage: %s <selected-subdiv.json> [level>=1] [out.obj] "
                     "[--camera-distance-scale s] [--fixed-edge-rate r]\n",
                     argv[0]);
        return 2;
    }

    const std::string inJson = argv[1];
    int level = 1;
    std::string outObj;
    float cameraDistanceScale = 1.0f;
    int fixedEdgeRate = 8;
    bool levelSet = false;
    bool outSet = false;
    for (int i = 2; i < argc; ++i)
    {
        const std::string arg = argv[i];
        if (arg == "--camera-distance-scale")
        {
            if (i + 1 >= argc)
            {
                std::fprintf(stderr, "Missing value for --camera-distance-scale\n");
                return 2;
            }
            cameraDistanceScale = std::strtof(argv[++i], nullptr);
            if (!(cameraDistanceScale > 0.0f))
            {
                std::fprintf(stderr, "camera-distance-scale must be > 0\n");
                return 2;
            }
            continue;
        }
        if (arg == "--fixed-edge-rate")
        {
            if (i + 1 >= argc)
            {
                std::fprintf(stderr, "Missing value for --fixed-edge-rate\n");
                return 2;
            }
            fixedEdgeRate = std::atoi(argv[++i]);
            if (fixedEdgeRate < 1)
            {
                std::fprintf(stderr, "fixed-edge-rate must be >= 1\n");
                return 2;
            }
            continue;
        }
        if (!levelSet)
        {
            level = std::max(1, std::atoi(argv[i]));
            levelSet = true;
            continue;
        }
        if (!outSet)
        {
            outObj = argv[i];
            outSet = true;
            continue;
        }
        std::fprintf(stderr,
                     "Usage: %s <selected-subdiv.json> [level>=1] [out.obj] "
                     "[--camera-distance-scale s] [--fixed-edge-rate r]\n",
                     argv[0]);
        return 2;
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
    const RefinedPositions refinedPositions = BuildRefinedPositions(*refiner, *patchTable, m.points);
    const std::vector<Primvar3> &patchEvalPositions = refinedPositions.values;
    const EdgeFactorCamera edgeFactorCamera =
        BuildEdgeFactorCamera(m.points, camera, cameraDistanceScale);
    const EdgeFactorSettings edgeFactorSettings = {};

    Far::PtexIndices ptex(*refiner);
    int ptexFaceCount = 0;
    std::vector<PtexFaceAdj> ptexFaceAdj;
    const int uniquePtexEdges =
        BuildUniquePtexEdgeIds(*refiner, level0, ptex, ptexFaceAdj, ptexFaceCount);
    int missingCenterPtexFaces = 0;
    const std::string missingNeighborsPath = outObj + ".missing_neighbors.txt";
    if (!WriteMissingPatchNeighbors(
            patchMap, ptexFaceAdj, ptexFaceCount, missingNeighborsPath, missingCenterPtexFaces))
    {
        std::fprintf(stderr, "Failed to write missing-patch neighbor report: %s\n", missingNeighborsPath.c_str());
        delete patchTable;
        delete refiner;
        return 1;
    }
    const EdgeFactorResult edgeFactors = ComputeEdgeFactors(patchMap,
                                                            *patchTable,
                                                            patchEvalPositions,
                                                            ptexFaceAdj,
                                                            ptexFaceCount,
                                                            uniquePtexEdges,
                                                            edgeFactorCamera,
                                                            edgeFactorSettings);

    int skippedPtexFaces = 0;
    const TriMesh fixedRateNoStitch = TessellateFixedRateNoStitch(
        patchMap, *patchTable, patchEvalPositions, ptexFaceCount, fixedEdgeRate, skippedPtexFaces);
    if (!WriteTriObj(outObj, fixedRateNoStitch))
    {
        std::fprintf(stderr, "Failed to write OBJ: %s\n", outObj.c_str());
        delete patchTable;
        delete refiner;
        return 1;
    }
    const std::string edgeRateObj = outObj + ".edge_rates.obj";
    const std::vector<ybi::testio::EdgeRateDebugLine> edgeRateLines =
        BuildEdgeRateDebugLines(
            patchMap, *patchTable, patchEvalPositions, ptexFaceAdj, ptexFaceCount, edgeFactors.edgeFactors);
    if (!ybi::testio::WriteEdgeRateDebugObj(
            edgeRateObj, edgeRateLines, edgeFactorSettings.minRate, edgeFactorSettings.maxRate))
    {
        std::fprintf(stderr, "Failed to write edge-rate debug OBJ: %s\n", edgeRateObj.c_str());
        delete patchTable;
        delete refiner;
        return 1;
    }

    std::printf("Wrote adaptive level-%d OBJ: %s\n", level, outObj.c_str());
    std::printf("Wrote edge-rate debug OBJ: %s\n", edgeRateObj.c_str());
    std::printf("Wrote missing-patch neighbors: %s\n", missingNeighborsPath.c_str());
    std::printf("  fixedRateNoStitch verts=%zu tris=%zu\n",
                fixedRateNoStitch.positions.size(),
                fixedRateNoStitch.indices.size() / 3);
    std::printf("  skippedPtexFaces=%d\n", skippedPtexFaces);
    std::printf("  missingCenterPtexFaces=%d\n", missingCenterPtexFaces);
    std::printf("  maxCalculatedEdgeFactor=%d\n", edgeFactors.maxCalculatedEdgeFactor);
    std::printf("  cameraDistanceScale=%g\n", cameraDistanceScale);
    std::printf("  fixedEdgeRate=%d\n", fixedEdgeRate);
    std::printf("  ptexFaces=%d\n", ptexFaceCount);
    std::printf("  controlCageEdgesWithOver2Faces=%d\n", edgesWithOver2Faces);

    delete patchTable;
    delete refiner;
    return 0;
}
