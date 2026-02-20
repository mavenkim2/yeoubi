#include "tesselation/diag_split_runtime.h"

#include "util/assert.h"

#include <opensubdiv/far/primvarRefiner.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace ybi::tesselation
{

namespace
{

pxr::GfVec3f Normalize(const pxr::GfVec3f &v)
{
    const float len = v.GetLength();
    if (len <= 1e-8f)
    {
        return pxr::GfVec3f(0.0f, 0.0f, 1.0f);
    }
    return v / len;
}

pxr::GfVec2f ToScreen(const EdgeFactorCamera &camera, const pxr::GfVec3f &point)
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

pxr::GfVec2f EdgeUV(int edge, float t)
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

} // namespace

Primvar3::Primvar3(const pxr::GfVec3f &p) : x(p[0]), y(p[1]), z(p[2]) {}

void Primvar3::Clear()
{
    x = y = z = 0.0f;
}

void Primvar3::AddWithWeight(const Primvar3 &p, float w)
{
    x += p.x * w;
    y += p.y * w;
    z += p.z * w;
}

pxr::GfVec3f Primvar3::ToGf() const
{
    return pxr::GfVec3f(x, y, z);
}

bool EvaluateLimitPosition(const OpenSubdiv::Far::PatchMap &patchMap,
                           const OpenSubdiv::Far::PatchTable &patchTable,
                           const std::vector<Primvar3> &positions,
                           int ptexFace,
                           float u,
                           float v,
                           pxr::GfVec3f &outPos)
{
    const OpenSubdiv::Far::PatchTable::PatchHandle *handle = patchMap.FindPatch(ptexFace, u, v);
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

    const OpenSubdiv::Far::ConstIndexArray cvIndices = patchTable.GetPatchVertices(*handle);
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

int ComputeDiagSplitEdgeFactor(const OpenSubdiv::Far::PatchMap &patchMap,
                               const OpenSubdiv::Far::PatchTable &patchTable,
                               const std::vector<Primvar3> &positions,
                               int ptexFace,
                               const pxr::GfVec2f &uvStart,
                               const pxr::GfVec2f &uvEnd,
                               const EdgeFactorCamera &camera,
                               const EdgeFactorSettings &settings)
{
    const int sampleCount = std::max(2, settings.numSamples);
    float maxLi = 0.0f;
    float sumLi = 0.0f;

    pxr::GfVec3f p0(0.0f);
    if (!EvaluateLimitPosition(
            patchMap, patchTable, positions, ptexFace, uvStart[0], uvStart[1], p0))
    {
        return settings.minRate;
    }
    pxr::GfVec2f prev = ToScreen(camera, p0);
    for (int i = 1; i < sampleCount; ++i)
    {
        const float t = float(i) / float(sampleCount - 1);
        const pxr::GfVec2f uv = uvStart * (1.0f - t) + uvEnd * t;
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

int ComputeDiagSplitEdgeFactor(const OpenSubdiv::Far::PatchMap &patchMap,
                               const OpenSubdiv::Far::PatchTable &patchTable,
                               const std::vector<Primvar3> &positions,
                               int ptexFace,
                               int edge,
                               const EdgeFactorCamera &camera,
                               const EdgeFactorSettings &settings)
{
    const pxr::GfVec2f uv0 = EdgeUV(edge, 0.0f);
    const pxr::GfVec2f uv1 = EdgeUV(edge, 1.0f);
    return ComputeDiagSplitEdgeFactor(
        patchMap, patchTable, positions, ptexFace, uv0, uv1, camera, settings);
}

EdgeFactorCamera BuildEdgeFactorCamera(const pxr::VtVec3fArray &points,
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

RefinedPositions BuildRefinedPositions(const OpenSubdiv::Far::TopologyRefiner &refiner,
                                       const OpenSubdiv::Far::PatchTable &patchTable,
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
    OpenSubdiv::Far::PrimvarRefiner prim(refiner);
    for (int l = 0; l < refiner.GetMaxLevel(); ++l)
    {
        Primvar3 *dst = src + refiner.GetLevel(l).GetNumVertices();
        prim.Interpolate(l + 1, src, dst);
        src = dst;
    }
    if (numLocalPoints > 0)
    {
        patchTable.ComputeLocalPointValues(out.values.data(),
                                           out.values.data() + out.numRefinerVerts);
    }
    return out;
}

int BuildUniquePtexEdgeIds(const OpenSubdiv::Far::TopologyRefiner &refiner,
                           const OpenSubdiv::Far::TopologyLevel &level0,
                           const OpenSubdiv::Far::PtexIndices &ptex,
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

EdgeFactorResult ComputeEdgeFactors(const OpenSubdiv::Far::PatchMap &patchMap,
                                    const OpenSubdiv::Far::PatchTable &patchTable,
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
            const int candidate = ComputeDiagSplitEdgeFactor(
                patchMap, patchTable, positions, pf, e, camera, settings);
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

void ApplyNgonNonUniformConstraint(const std::vector<PtexFaceAdj> &faces,
                                   int numPtexFaces,
                                   std::vector<int> &edgeRates)
{
    for (int pf = 0; pf < numPtexFaces; ++pf)
    {
        for (int e = 0; e < 4; ++e)
        {
            const int af = faces[pf].adjFace[e];
            if (af < 0 || af >= numPtexFaces)
            {
                continue;
            }
            const int edgeId = faces[pf].edgeIndex[e];
            if (edgeId < 0 || edgeId >= int(edgeRates.size()))
            {
                continue;
            }
            if (faces[pf].fromNgon || faces[af].fromNgon)
            {
                edgeRates[edgeId] = kDiagSplitNonUniform;
            }
        }
    }
}

std::vector<ybi::testio::EdgeRateDebugLine>
BuildEdgeRateDebugLines(const OpenSubdiv::Far::PatchMap &patchMap,
                        const OpenSubdiv::Far::PatchTable &patchTable,
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
                if (!EvaluateLimitPosition(patchMap, patchTable, positions, pf, uv[0], uv[1], p))
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

} // namespace ybi::tesselation
