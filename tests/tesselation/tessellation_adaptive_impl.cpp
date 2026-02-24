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

static bool EvaluateLimitPosition(const Far::PatchMap &patchMap,
                                  const Far::PatchTable &patchTable,
                                  const std::vector<LimitEvalVertex> &limitValues,
                                  int ptexFaceId,
                                  const pxr::GfVec2f &uv,
                                  pxr::GfVec3f *outP)
{
    if (!outP)
    {
        return false;
    }
    const Far::PatchTable::PatchHandle *handle = patchMap.FindPatch(ptexFaceId, uv[0], uv[1]);
    if (!handle)
    {
        return false;
    }

    float pWeights[20] = {0.0f};
    patchTable.EvaluateBasis(*handle, uv[0], uv[1], pWeights);
    Far::ConstIndexArray cvs = patchTable.GetPatchVertices(*handle);

    pxr::GfVec3f p(0.0f);
    for (int i = 0; i < cvs.size(); ++i)
    {
        p += limitValues[cvs[i]].p * pWeights[i];
    }
    *outP = p;
    return true;
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

static int ComputeDiagSplitPatchEdgeFactor(const Far::PatchMap &patchMap,
                                           const Far::PatchTable &patchTable,
                                           const std::vector<LimitEvalVertex> &limitValues,
                                           int ptexFaceId,
                                           const pxr::GfVec2f &uvStart,
                                           const pxr::GfVec2f &uvEnd,
                                           int sampleStepsN,
                                           float targetPixelSpacing,
                                           int splitThreshold,
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
    float sumLi = 0.0f;
    pxr::GfVec3f p0(0.0f);
    if (!EvaluateLimitPosition(patchMap, patchTable, limitValues, ptexFaceId, uvStart, &p0))
    {
        return SUBDIV_EDGE_FACTOR_NON_UNIFORM;
    }
    pxr::GfVec2f prev =
        ProjectToScreen(p0, eye, lookAt, viewportWidth, viewportHeight, verticalFovDegrees);
    for (int i = 1; i < sampleStepsN; ++i)
    {
        const float t = float(i) / float(sampleStepsN - 1);
        const pxr::GfVec2f uv = uvStart * (1.0f - t) + uvEnd * t;
        pxr::GfVec3f p(0.0f);
        if (!EvaluateLimitPosition(patchMap, patchTable, limitValues, ptexFaceId, uv, &p))
        {
            return SUBDIV_EDGE_FACTOR_NON_UNIFORM;
        }
        const pxr::GfVec2f s =
            ProjectToScreen(p, eye, lookAt, viewportWidth, viewportHeight, verticalFovDegrees);
        const pxr::GfVec2f d = s - prev;
        const float li = std::sqrt(d[0] * d[0] + d[1] * d[1]);
        sumLi += li;
        maxLi = std::max(maxLi, li);
        prev = s;
    }

    const int tMin = std::max(1, int(std::ceil(sumLi / targetPixelSpacing)));
    const int tMax =
        std::max(1, int(std::ceil(float(sampleStepsN - 1) * maxLi / targetPixelSpacing)));
    if ((tMax - tMin) >= splitThreshold)
    {
        return SUBDIV_EDGE_FACTOR_NON_UNIFORM;
    }
    return tMax;
}

static SubdivisionEdge &GetOrCreateEdge(SubdivisionEdgeMap &edgeMap,
                                        int v0,
                                        int v1,
                                        int ptexFaceId,
                                        const pxr::GfVec2f &uv0,
                                        const pxr::GfVec2f &uv1)
{
    const uint64_t key = MakeEdgeKey(v0, v1);
    auto it = edgeMap.find(key);
    if (it == edgeMap.end())
    {
        SubdivisionEdge edge = {};
        edge.v0 = v0;
        edge.v1 = v1;
        edge.sampleVStart = v0;
        edge.sampleVEnd = v1;
        edge.storedPtexFaceId = ptexFaceId;
        edge.storedUv0 = uv0;
        edge.storedUv1 = uv1;
        edge.hasStoredPatchParams = (ptexFaceId >= 0);
        it = edgeMap.emplace(key, edge).first;
    }
    return it->second;
}

static void SetEdgeSampleParams(SubdivisionEdge &edge,
                                int vStart,
                                int vEnd,
                                int ptexFaceId,
                                const pxr::GfVec2f &uvStart,
                                const pxr::GfVec2f &uvEnd,
                                bool overwrite)
{
    if (!overwrite && edge.hasStoredPatchParams && edge.sampleVStart >= 0 && edge.sampleVEnd >= 0)
    {
        return;
    }
    edge.sampleVStart = vStart;
    edge.sampleVEnd = vEnd;
    edge.storedPtexFaceId = ptexFaceId;
    edge.storedUv0 = uvStart;
    edge.storedUv1 = uvEnd;
    edge.hasStoredPatchParams = (ptexFaceId >= 0);
}

static int GetEdgeFactor(const SubdivisionEdgeMap &edgeMap, int v0, int v1)
{
    const uint64_t key = MakeEdgeKey(v0, v1);
    const auto it = edgeMap.find(key);
    YBI_ERROR(it != edgeMap.end(), "Missing edge factor lookup for edge (%d, %d)\n", v0, v1);
    const int t = it->second.tmaxEdgeFactor;
    return (t >= 1) ? std::min(t, 2) : t;
}

static int EnsurePatchEdgeFactor(const SelectedSubdivMesh &m,
                                 const SubdivisionPatch &patch,
                                 int edgeIndex,
                                 SubdivisionEdgeMap &edgeMap,
                                 int &nextGeneratedVertexId,
                                 bool allowQuadCoarseUpgrade,
                                 const Far::PatchMap &patchMap,
                                 const Far::PatchTable &patchTable,
                                 const std::vector<LimitEvalVertex> &limitValues,
                                 int sampleStepsN,
                                 float targetPixelSpacing,
                                 int splitThreshold,
                                 const pxr::GfVec3f &eye,
                                 const pxr::GfVec3f &lookAt,
                                 int viewportWidth,
                                 int viewportHeight,
                                 float verticalFovDegrees,
                                 int *computedCountOut)
{
    const int next = (edgeIndex + 1) & 3;
    SubdivisionEdge &edge = GetOrCreateEdge(edgeMap,
                                            patch.verts[edgeIndex],
                                            patch.verts[next],
                                            patch.ptexFaceId,
                                            patch.uv[edgeIndex],
                                            patch.uv[next]);
    if (edge.tmaxEdgeFactor == SUBDIV_EDGE_FACTOR_UNINITIALIZED)
    {
        edge.tmaxEdgeFactor = ComputeDiagSplitPatchEdgeFactor(patchMap,
                                                              patchTable,
                                                              limitValues,
                                                              patch.ptexFaceId,
                                                              patch.uv[edgeIndex],
                                                              patch.uv[next],
                                                              sampleStepsN,
                                                              targetPixelSpacing,
                                                              splitThreshold,
                                                              eye,
                                                              lookAt,
                                                              viewportWidth,
                                                              viewportHeight,
                                                              verticalFovDegrees);
        if (edge.tmaxEdgeFactor != SUBDIV_EDGE_FACTOR_NON_UNIFORM && edge.tmaxEdgeFactor < 1)
        {
            edge.tmaxEdgeFactor = 1;
        }
        if (edge.tmaxEdgeFactor >= 1)
        {
            edge.tmaxEdgeFactor = std::min(edge.tmaxEdgeFactor, 2);
        }
        const bool coarseFaceIsQuad = (patch.coarseFace >= 0) &&
                                      (patch.coarseFace < int(m.faceVertexCounts.size())) &&
                                      (m.faceVertexCounts[patch.coarseFace] == 4);
        if (allowQuadCoarseUpgrade && coarseFaceIsQuad && edge.tmaxEdgeFactor == 1)
        {
            edge.tmaxEdgeFactor = 2;
        }
        if (edge.tmaxEdgeFactor >= 1)
        {
            edge.transitionedUninitializedToUniform = true;
            SetEdgeSampleParams(edge,
                                patch.verts[edgeIndex],
                                patch.verts[next],
                                patch.ptexFaceId,
                                patch.uv[edgeIndex],
                                patch.uv[next],
                                /*overwrite*/ true);
            if (edge.tmaxEdgeFactor >= 2)
            {
                const int splitPosStored = edge.tmaxEdgeFactor / 2;
                edge.edgeVertexIndexStart = nextGeneratedVertexId;
                nextGeneratedVertexId += (edge.tmaxEdgeFactor - 1);
                edge.midpointVertex = edge.edgeVertexIndexStart + (splitPosStored - 1);
            }
            else
            {
                edge.edgeVertexIndexStart = -1;
                edge.midpointVertex = -1;
            }
        }
        else
        {
            edge.edgeVertexIndexStart = -1;
        }
        if (computedCountOut)
        {
            (*computedCountOut)++;
        }
    }
    else if (edge.tmaxEdgeFactor >= 1)
    {
        edge.tmaxEdgeFactor = std::min(edge.tmaxEdgeFactor, 2);
        SetEdgeSampleParams(edge,
                            patch.verts[edgeIndex],
                            patch.verts[next],
                            patch.ptexFaceId,
                            patch.uv[edgeIndex],
                            patch.uv[next],
                            /*overwrite*/ false);
    }
    return edge.tmaxEdgeFactor;
}

static int PrecomputePatchEdgeFactors(const SelectedSubdivMesh &m,
                                      const std::vector<SubdivisionPatch> &patches,
                                      SubdivisionEdgeMap &edgeMap,
                                      int &nextGeneratedVertexId,
                                      const Far::PatchMap &patchMap,
                                      const Far::PatchTable &patchTable,
                                      const std::vector<LimitEvalVertex> &limitValues,
                                      int sampleStepsN,
                                      float targetPixelSpacing,
                                      int splitThreshold,
                                      const pxr::GfVec3f &eye,
                                      const pxr::GfVec3f &lookAt,
                                      int viewportWidth,
                                      int viewportHeight,
                                      float verticalFovDegrees)
{
    int computedCount = 0;
    for (const SubdivisionPatch &patch : patches)
    {
        for (int e = 0; e < 4; ++e)
        {
            EnsurePatchEdgeFactor(m,
                                  patch,
                                  e,
                                  edgeMap,
                                  nextGeneratedVertexId,
                                  /*allowQuadCoarseUpgrade*/ true,
                                  patchMap,
                                  patchTable,
                                  limitValues,
                                  sampleStepsN,
                                  targetPixelSpacing,
                                  splitThreshold,
                                  eye,
                                  lookAt,
                                  viewportWidth,
                                  viewportHeight,
                                  verticalFovDegrees,
                                  &computedCount);
        }
    }
    return computedCount;
}

static std::vector<SubdivisionPatch>
DiagSplitPatches(const SelectedSubdivMesh &m,
                 const std::vector<SubdivisionPatch> &patches,
                 SubdivisionEdgeMap &edgeMap,
                 int &nextGeneratedVertexId,
                 const Far::PatchMap &patchMap,
                 const Far::PatchTable &patchTable,
                 const std::vector<LimitEvalVertex> &limitValues,
                 int sampleStepsN,
                 float targetPixelSpacing,
                 int splitThreshold,
                 const pxr::GfVec3f &eye,
                 const pxr::GfVec3f &lookAt,
                 int viewportWidth,
                 int viewportHeight,
                 float verticalFovDegrees,
                 int *computedCountOut)
{
    std::vector<SubdivisionPatch> worklist = patches;
    std::vector<SubdivisionPatch> out;
    out.reserve(patches.size());

    size_t patchCursor = 0;
    while (patchCursor < worklist.size())
    {
        const SubdivisionPatch patch = worklist[patchCursor++];

        int edgeFactor[4] = {SUBDIV_EDGE_FACTOR_UNINITIALIZED,
                             SUBDIV_EDGE_FACTOR_UNINITIALIZED,
                             SUBDIV_EDGE_FACTOR_UNINITIALIZED,
                             SUBDIV_EDGE_FACTOR_UNINITIALIZED};
        int splitEdge = -1;
        for (int e = 0; e < 4; ++e)
        {
            edgeFactor[e] = EnsurePatchEdgeFactor(m,
                                                  patch,
                                                  e,
                                                  edgeMap,
                                                  nextGeneratedVertexId,
                                                  /*allowQuadCoarseUpgrade*/ false,
                                                  patchMap,
                                                  patchTable,
                                                  limitValues,
                                                  sampleStepsN,
                                                  targetPixelSpacing,
                                                  splitThreshold,
                                                  eye,
                                                  lookAt,
                                                  viewportWidth,
                                                  viewportHeight,
                                                  verticalFovDegrees,
                                                  computedCountOut);
            if (splitEdge < 0 && edgeFactor[e] == SUBDIV_EDGE_FACTOR_NON_UNIFORM)
            {
                splitEdge = e;
            }
        }

        if (splitEdge < 0)
        {
            out.push_back(patch);
            continue;
        }

        const int oppositeEdge = (splitEdge + 2) & 3;
        if (edgeFactor[oppositeEdge] == 1)
        {
            const int ov0 = patch.verts[oppositeEdge];
            const int ov1 = patch.verts[(oppositeEdge + 1) & 3];
            const int ovn = (oppositeEdge + 1) & 3;
            SubdivisionEdge &opp = GetOrCreateEdge(
                edgeMap, ov0, ov1, patch.ptexFaceId, patch.uv[oppositeEdge], patch.uv[ovn]);
            opp.tmaxEdgeFactor = 2;
            SetEdgeSampleParams(
                opp, ov0, ov1, patch.ptexFaceId, patch.uv[oppositeEdge], patch.uv[ovn], true);
            if (opp.edgeVertexIndexStart < 0 || opp.midpointVertex < 0)
            {
                opp.edgeVertexIndexStart = nextGeneratedVertexId++;
                opp.midpointVertex = opp.edgeVertexIndexStart;
            }
            edgeFactor[oppositeEdge] = 2;
        }

        const int splitPair[2] = {splitEdge, oppositeEdge};
        int splitVerts[2] = {-1, -1};
        pxr::GfVec2f splitUV[2];
        for (int i = 0; i < 2; ++i)
        {
            const int e = splitPair[i];
            const int n = (e + 1) & 3;
            const int v0 = patch.verts[e];
            const int v1 = patch.verts[n];
            SubdivisionEdge &edge =
                GetOrCreateEdge(edgeMap, v0, v1, patch.ptexFaceId, patch.uv[e], patch.uv[n]);

            int mid = edge.midpointVertex;
            int t = edgeFactor[e];
            float alpha = 0.5f;

            if (t != SUBDIV_EDGE_FACTOR_NON_UNIFORM)
            {
                YBI_ASSERT(t > 1);
                const int splitPosStored = t / 2;
                alpha = float(splitPosStored) / float(t);
                splitUV[i] = patch.uv[e] * (1.0f - alpha) + patch.uv[n] * alpha;
                YBI_ASSERT(mid >= 0);
                YBI_ASSERT(edge.edgeVertexIndexStart != -1);

                if (!edge.split)
                {
                    edge.split = true;

                    // NOTE: these are technically wrong when child edge factor is 1, but the
                    // lookup won't use edgeVertexIndexStart in that case.
                    SubdivisionEdge &edgeA = GetOrCreateEdge(
                        edgeMap, v0, mid, patch.ptexFaceId, patch.uv[e], splitUV[i]);
                    SubdivisionEdge &edgeB = GetOrCreateEdge(
                        edgeMap, mid, v1, patch.ptexFaceId, splitUV[i], patch.uv[n]);
                    edgeA.tmaxEdgeFactor = splitPosStored;
                    edgeA.edgeVertexIndexStart = edge.edgeVertexIndexStart;
                    edgeA.midpointVertex =
                        edgeA.tmaxEdgeFactor > 1
                            ? edgeA.edgeVertexIndexStart + (edgeA.tmaxEdgeFactor / 2) - 1
                            : -1;
                    SetEdgeSampleParams(
                        edgeA, v0, mid, patch.ptexFaceId, patch.uv[e], splitUV[i], true);

                    edgeB.tmaxEdgeFactor = t - splitPosStored;
                    edgeB.edgeVertexIndexStart = edge.edgeVertexIndexStart + splitPosStored;
                    edgeB.midpointVertex =
                        edgeB.tmaxEdgeFactor > 1
                            ? edgeB.edgeVertexIndexStart + (edgeB.tmaxEdgeFactor / 2) - 1
                            : -1;
                    SetEdgeSampleParams(
                        edgeB, mid, v1, patch.ptexFaceId, splitUV[i], patch.uv[n], true);

                    edgeA.boundary = edge.boundary;
                    edgeB.boundary = edge.boundary;
                }
            }
            else
            {
                splitUV[i] = patch.uv[e] * (1.0f - alpha) + patch.uv[n] * alpha;
                if (!edge.split)
                {
                    edge.split = true;
                    mid = nextGeneratedVertexId++;
                    edge.midpointVertex = mid;
                    SubdivisionEdge &edgeA = GetOrCreateEdge(
                        edgeMap, v0, mid, patch.ptexFaceId, patch.uv[e], splitUV[i]);
                    SubdivisionEdge &edgeB = GetOrCreateEdge(
                        edgeMap, mid, v1, patch.ptexFaceId, splitUV[i], patch.uv[n]);
                    edgeA.tmaxEdgeFactor = SUBDIV_EDGE_FACTOR_UNINITIALIZED;
                    edgeB.tmaxEdgeFactor = SUBDIV_EDGE_FACTOR_UNINITIALIZED;
                    edgeA.boundary = edge.boundary;
                    edgeB.boundary = edge.boundary;
                }
            }

            splitVerts[i] = mid;
        }

        SubdivisionEdge &seam = GetOrCreateEdge(
            edgeMap, splitVerts[0], splitVerts[1], patch.ptexFaceId, splitUV[0], splitUV[1]);
        seam.tmaxEdgeFactor = SUBDIV_EDGE_FACTOR_UNINITIALIZED;
        seam.edgeVertexIndexStart = -1;
        seam.boundary = false;

        const int i0 = splitEdge;
        const int i1 = (i0 + 1) & 3;
        const int i2 = (i0 + 2) & 3;
        const int i3 = (i0 + 3) & 3;
        const int midA = splitVerts[0];
        const int midB = splitVerts[1];

        SubdivisionPatch childA = patch;
        childA.verts[0] = patch.verts[i0];
        childA.verts[1] = midA;
        childA.verts[2] = midB;
        childA.verts[3] = patch.verts[i3];
        childA.uv[0] = patch.uv[i0];
        childA.uv[1] = splitUV[0];
        childA.uv[2] = splitUV[1];
        childA.uv[3] = patch.uv[i3];

        SubdivisionPatch childB = patch;
        childB.verts[0] = midA;
        childB.verts[1] = patch.verts[i1];
        childB.verts[2] = patch.verts[i2];
        childB.verts[3] = midB;
        childB.uv[0] = splitUV[0];
        childB.uv[1] = patch.uv[i1];
        childB.uv[2] = patch.uv[i2];
        childB.uv[3] = splitUV[1];

        worklist.push_back(childA);
        worklist.push_back(childB);
    }

    return out;
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
            pxr::GfVec2f uv0(0.0f, 0.0f);
            pxr::GfVec2f uv1(0.0f, 0.0f);
            if (n == 4)
            {
                const pxr::GfVec2f cornerUVs[4] = {pxr::GfVec2f(0.0f, 0.0f),
                                                   pxr::GfVec2f(1.0f, 0.0f),
                                                   pxr::GfVec2f(1.0f, 1.0f),
                                                   pxr::GfVec2f(0.0f, 1.0f)};
                uv0 = cornerUVs[i];
                uv1 = cornerUVs[(i + 1) % n];
            }
            else
            {
                uv0 = pxr::GfVec2f(float(i) / float(n), 0.0f);
                uv1 = pxr::GfVec2f(float(i + 1) / float(n), 0.0f);
            }
            SubdivisionEdge &edge = GetOrCreateEdge(edgeMap, a, b, int(f), uv0, uv1);
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

static int GetOrAllocateMidpointVertex(SubdivisionEdgeMap &edgeMap,
                                       int v0,
                                       int v1,
                                       int ptexFaceId,
                                       const pxr::GfVec2f &uv0,
                                       const pxr::GfVec2f &uv1,
                                       int &nextVertexId)
{
    SubdivisionEdge &edge = GetOrCreateEdge(edgeMap, v0, v1, ptexFaceId, uv0, uv1);
    if (edge.midpointVertex < 0)
    {
        edge.midpointVertex = nextVertexId++;
    }
    return edge.midpointVertex;
}

static void AddEdgeUse(SubdivisionEdgeMap &edgeMap,
                       int v0,
                       int v1,
                       int faceId,
                       int ptexFaceId,
                       const pxr::GfVec2f &uv0,
                       const pxr::GfVec2f &uv1)
{
    SubdivisionEdge &edge = GetOrCreateEdge(edgeMap, v0, v1, ptexFaceId, uv0, uv1);
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
            patch.uv[0] = pxr::GfVec2f(0.0f, 0.0f);
            patch.uv[1] = pxr::GfVec2f(1.0f, 0.0f);
            patch.uv[2] = pxr::GfVec2f(1.0f, 1.0f);
            patch.uv[3] = pxr::GfVec2f(0.0f, 1.0f);
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
                const int midNext = GetOrAllocateMidpointVertex(edgeMap,
                                                                v,
                                                                vNext,
                                                                basePtexFaceId + i,
                                                                pxr::GfVec2f(0.0f, 0.0f),
                                                                pxr::GfVec2f(1.0f, 0.0f),
                                                                nextVertexId);
                const int midPrev = GetOrAllocateMidpointVertex(edgeMap,
                                                                vPrev,
                                                                v,
                                                                basePtexFaceId + i,
                                                                pxr::GfVec2f(0.0f, 1.0f),
                                                                pxr::GfVec2f(0.0f, 0.0f),
                                                                nextVertexId);

                SubdivisionPatch patch = {};
                patch.verts[0] = v;
                patch.verts[1] = midNext;
                patch.verts[2] = center;
                patch.verts[3] = midPrev;
                patch.uv[0] = pxr::GfVec2f(0.0f, 0.0f);
                patch.uv[1] = pxr::GfVec2f(1.0f, 0.0f);
                patch.uv[2] = pxr::GfVec2f(1.0f, 1.0f);
                patch.uv[3] = pxr::GfVec2f(0.0f, 1.0f);
                patch.coarseFace = int(f);
                patch.quadrant = i;
                patch.ptexFaceId = basePtexFaceId + i;
                const int patchId = int(patches.size());
                patches.push_back(patch);

                // Only generated edges for non-quad quadrangulation.
                AddEdgeUse(
                    edgeMap, v, midNext, patchId, patch.ptexFaceId, patch.uv[0], patch.uv[1]);
                AddEdgeUse(
                    edgeMap, midNext, center, patchId, patch.ptexFaceId, patch.uv[1], patch.uv[2]);
                AddEdgeUse(
                    edgeMap, center, midPrev, patchId, patch.ptexFaceId, patch.uv[2], patch.uv[3]);
                AddEdgeUse(
                    edgeMap, midPrev, v, patchId, patch.ptexFaceId, patch.uv[3], patch.uv[0]);
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
        if (it.second.tmaxEdgeFactor != SUBDIV_EDGE_FACTOR_UNINITIALIZED)
        {
            count++;
        }
    }
    return count;
}

static bool VerifyInitializedUniformEdgesHaveStoredPatchParams(const SubdivisionEdgeMap &edgeMap,
                                                               int *missingStoredParamsOut,
                                                               int *badUniformFactorOut)
{
    int missingStoredParams = 0;
    int badUniformFactor = 0;
    for (const auto &it : edgeMap)
    {
        const SubdivisionEdge &edge = it.second;
        if (!edge.transitionedUninitializedToUniform)
        {
            continue;
        }
        if (edge.tmaxEdgeFactor < 1)
        {
            badUniformFactor++;
        }
        if (!edge.hasStoredPatchParams || edge.storedPtexFaceId < 0 || edge.sampleVStart < 0 ||
            edge.sampleVEnd < 0)
        {
            missingStoredParams++;
        }
    }

    if (missingStoredParamsOut)
    {
        *missingStoredParamsOut = missingStoredParams;
    }
    if (badUniformFactorOut)
    {
        *badUniformFactorOut = badUniformFactor;
    }
    return (missingStoredParams == 0) && (badUniformFactor == 0);
}

static bool WriteLeafPatchInnerGridObj(const std::vector<SubdivisionPatch> &leafPatches,
                                       const SubdivisionEdgeMap &edgeMap,
                                       const Far::PatchMap &patchMap,
                                       const Far::PatchTable &patchTable,
                                       const std::vector<LimitEvalVertex> &limitValues,
                                       const std::string &outObjPath,
                                       int *outVertexCount,
                                       int *outTriangleCount)
{
    if (outVertexCount)
    {
        *outVertexCount = 0;
    }
    if (outTriangleCount)
    {
        *outTriangleCount = 0;
    }

    FILE *f = std::fopen(outObjPath.c_str(), "w");
    if (!f)
    {
        return false;
    }

    int maxVertexId = -1;
    for (const SubdivisionPatch &patch : leafPatches)
    {
        for (int i = 0; i < 4; ++i)
        {
            maxVertexId = std::max(maxVertexId, patch.verts[i]);
        }
    }
    for (const auto &it : edgeMap)
    {
        const SubdivisionEdge &edge = it.second;
        maxVertexId = std::max(maxVertexId, edge.v0);
        maxVertexId = std::max(maxVertexId, edge.v1);
        maxVertexId = std::max(maxVertexId, edge.sampleVStart);
        maxVertexId = std::max(maxVertexId, edge.sampleVEnd);
        maxVertexId = std::max(maxVertexId, edge.midpointVertex);
        if (edge.edgeVertexIndexStart >= 0 && edge.tmaxEdgeFactor >= 2)
        {
            maxVertexId = std::max(
                maxVertexId, edge.edgeVertexIndexStart + (std::min(edge.tmaxEdgeFactor, 2) - 2));
        }
    }
    if (maxVertexId < 0)
    {
        std::fclose(f);
        return true;
    }

    const float inf = std::numeric_limits<float>::infinity();
    std::vector<pxr::GfVec3f> sharedPos(maxVertexId + 1, pxr::GfVec3f(inf, inf, inf));
    std::vector<uint8_t> sharedInit(maxVertexId + 1, uint8_t(0));
    auto isSharedInit = [&](int id) -> bool {
        return (id >= 0) && (id < int(sharedInit.size())) && (sharedInit[id] != 0);
    };

    for (const auto &it : edgeMap)
    {
        const SubdivisionEdge &edge = it.second;
        if (!edge.transitionedUninitializedToUniform)
        {
            continue;
        }

        YBI_ASSERT(edge.hasStoredPatchParams);
        YBI_ASSERT(edge.storedPtexFaceId >= 0);
        YBI_ASSERT(edge.sampleVStart >= 0 && edge.sampleVStart < int(sharedPos.size()));
        YBI_ASSERT(edge.sampleVEnd >= 0 && edge.sampleVEnd < int(sharedPos.size()));

        int t = edge.tmaxEdgeFactor;
        YBI_ASSERT(t >= 1);
        t = std::min(t, 2);
        if (t > 1)
        {
            YBI_ASSERT(edge.edgeVertexIndexStart >= 0);
            YBI_ASSERT(edge.edgeVertexIndexStart < int(sharedPos.size()));
        }

        for (int k = 0; k <= t; ++k)
        {
            const float a = float(k) / float(t);
            const pxr::GfVec2f uv = edge.storedUv0 * (1.0f - a) + edge.storedUv1 * a;
            pxr::GfVec3f p(0.0f);
            if (!EvaluateLimitPosition(
                    patchMap, patchTable, limitValues, edge.storedPtexFaceId, uv, &p))
            {
                std::fclose(f);
                return false;
            }

            int vertexId = -1;
            if (k == 0)
            {
                vertexId = edge.sampleVStart;
            }
            else if (k == t)
            {
                vertexId = edge.sampleVEnd;
            }
            else
            {
                vertexId = edge.edgeVertexIndexStart + (k - 1);
            }
            YBI_ASSERT(vertexId >= 0 && vertexId < int(sharedPos.size()));
            if (!isSharedInit(vertexId))
            {
                sharedPos[vertexId] = p;
                sharedInit[vertexId] = 1;
            }
        }
    }

    int vertexCount = 0;
    int triangleCount = 0;
    std::vector<int> sharedObjIndex(maxVertexId + 1, 0);
    for (int i = 0; i < int(sharedPos.size()); ++i)
    {
        if (!isSharedInit(i))
        {
            continue;
        }
        const pxr::GfVec3f &p = sharedPos[i];
        std::fprintf(f, "v %.9g %.9g %.9g\n", double(p[0]), double(p[1]), double(p[2]));
        vertexCount++;
        sharedObjIndex[i] = vertexCount;
    }

    auto emitTri = [&](int a, int b, int c) {
        YBI_ASSERT(a > 0 && b > 0 && c > 0);
        if (a == b || b == c || a == c)
        {
            return;
        }
        std::fprintf(f, "f %d %d %d\n", a, b, c);
        triangleCount++;
    };

    auto getEdgeRef = [&](int a, int b) -> const SubdivisionEdge & {
        const uint64_t key = MakeEdgeKey(a, b);
        const auto it = edgeMap.find(key);
        YBI_ERROR(it != edgeMap.end(), "Missing edge for stitching (%d, %d)\n", a, b);
        return it->second;
    };

    auto getOrientedEdgeVertexId =
        [&](const SubdivisionEdge &edge, int a, int b, int t, int k) -> int {
        YBI_ASSERT(t >= 1);
        YBI_ASSERT(k >= 0 && k <= t);
        const bool forward = (edge.v0 == a && edge.v1 == b);
        const bool reverse = (edge.v0 == b && edge.v1 == a);
        YBI_ASSERT(forward || reverse);
        if (k == 0)
        {
            return a;
        }
        if (k == t)
        {
            return b;
        }
        YBI_ASSERT(t > 1);
        YBI_ASSERT(edge.edgeVertexIndexStart >= 0);
        const int localOffset = forward ? (k - 1) : (t - 1 - k);
        return edge.edgeVertexIndexStart + localOffset;
    };

    auto stitchStrip = [&](const std::vector<int> &outer, const std::vector<int> &inner) {
        YBI_ASSERT(!outer.empty());
        YBI_ASSERT(!inner.empty());

        size_t io = 0;
        size_t ii = 0;
        if (inner.size() == 1)
        {
            while (io + 1 < outer.size())
            {
                emitTri(outer[io], outer[io + 1], inner[0]);
                io++;
            }
            return;
        }

        while (io + 1 < outer.size() || ii + 1 < inner.size())
        {
            const bool canAdvanceOuter = (io + 1 < outer.size());
            const bool canAdvanceInner = (ii + 1 < inner.size());

            if (!canAdvanceInner)
            {
                emitTri(outer[io], outer[io + 1], inner[ii]);
                io++;
                continue;
            }
            if (!canAdvanceOuter)
            {
                emitTri(outer[io], inner[ii + 1], inner[ii]);
                ii++;
                continue;
            }

            const float uOuter = float(io + 1) / float(outer.size() - 1);
            const float uInner = float(ii + 1) / float(inner.size() - 1);
            if (uOuter <= uInner)
            {
                emitTri(outer[io], outer[io + 1], inner[ii]);
                io++;
            }
            else
            {
                emitTri(outer[io], inner[ii + 1], inner[ii]);
                ii++;
            }
        }
    };

    for (const SubdivisionPatch &patch : leafPatches)
    {
        const int e0 = GetEdgeFactor(edgeMap, patch.verts[0], patch.verts[1]);
        const int e1 = GetEdgeFactor(edgeMap, patch.verts[1], patch.verts[2]);
        const int e2 = GetEdgeFactor(edgeMap, patch.verts[2], patch.verts[3]);
        const int e3 = GetEdgeFactor(edgeMap, patch.verts[3], patch.verts[0]);

        YBI_ASSERT(e0 != SUBDIV_EDGE_FACTOR_NON_UNIFORM);
        YBI_ASSERT(e1 != SUBDIV_EDGE_FACTOR_NON_UNIFORM);
        YBI_ASSERT(e2 != SUBDIV_EDGE_FACTOR_NON_UNIFORM);
        YBI_ASSERT(e3 != SUBDIV_EDGE_FACTOR_NON_UNIFORM);
        YBI_ASSERT(e0 >= 1 && e1 >= 1 && e2 >= 1 && e3 >= 1);

        const int nu = std::max(std::max(e0, e2), 2);
        const int nv = std::max(std::max(e1, e3), 2);
        const int cols = nu - 1;
        const int rows = nv - 1;
        std::vector<int> innerGridObj(cols * rows, 0);
        auto innerAt = [&](int iu, int iv) -> int & {
            YBI_ASSERT(iu >= 1 && iu <= cols);
            YBI_ASSERT(iv >= 1 && iv <= rows);
            return innerGridObj[(iv - 1) * cols + (iu - 1)];
        };

        for (int iv = 1; iv < nv; ++iv)
        {
            const float sv = float(iv) / float(nv);
            for (int iu = 1; iu < nu; ++iu)
            {
                const float su = float(iu) / float(nu);
                const pxr::GfVec2f uvBottom = patch.uv[0] * (1.0f - su) + patch.uv[1] * su;
                const pxr::GfVec2f uvTop = patch.uv[3] * (1.0f - su) + patch.uv[2] * su;
                const pxr::GfVec2f uv = uvBottom * (1.0f - sv) + uvTop * sv;

                pxr::GfVec3f p(0.0f);
                if (!EvaluateLimitPosition(
                        patchMap, patchTable, limitValues, patch.ptexFaceId, uv, &p))
                {
                    std::fclose(f);
                    return false;
                }
                std::fprintf(f, "v %.9g %.9g %.9g\n", double(p[0]), double(p[1]), double(p[2]));
                vertexCount++;
                innerAt(iu, iv) = vertexCount;
            }
        }

        if (cols >= 2 && rows >= 2)
        {
            for (int y = 0; y < rows - 1; ++y)
            {
                for (int x = 0; x < cols - 1; ++x)
                {
                    const int i00 = innerAt(x + 1, y + 1);
                    const int i10 = innerAt(x + 2, y + 1);
                    const int i01 = innerAt(x + 1, y + 2);
                    const int i11 = innerAt(x + 2, y + 2);
                    emitTri(i00, i10, i11);
                    emitTri(i00, i11, i01);
                }
            }
        }

        const SubdivisionEdge &edge0 = getEdgeRef(patch.verts[0], patch.verts[1]);
        const SubdivisionEdge &edge1 = getEdgeRef(patch.verts[1], patch.verts[2]);
        const SubdivisionEdge &edge2 = getEdgeRef(patch.verts[2], patch.verts[3]);
        const SubdivisionEdge &edge3 = getEdgeRef(patch.verts[3], patch.verts[0]);

        auto buildOuter = [&](const SubdivisionEdge &edge, int a, int b, int t) {
            std::vector<int> outer;
            outer.reserve(t + 1);
            for (int k = 0; k <= t; ++k)
            {
                const int vertexId = getOrientedEdgeVertexId(edge, a, b, t, k);
                YBI_ASSERT(vertexId >= 0 && vertexId < int(sharedPos.size()));
                YBI_ASSERT(isSharedInit(vertexId));
                const int objId = sharedObjIndex[vertexId];
                YBI_ASSERT(objId > 0);
                outer.push_back(objId);
            }
            return outer;
        };

        std::vector<int> outer0 = buildOuter(edge0, patch.verts[0], patch.verts[1], e0);
        std::vector<int> outer1 = buildOuter(edge1, patch.verts[1], patch.verts[2], e1);
        std::vector<int> outer2 = buildOuter(edge2, patch.verts[2], patch.verts[3], e2);
        std::vector<int> outer3 = buildOuter(edge3, patch.verts[3], patch.verts[0], e3);

        std::vector<int> inner0;
        std::vector<int> inner1;
        std::vector<int> inner2;
        std::vector<int> inner3;
        inner0.reserve(cols);
        inner1.reserve(rows);
        inner2.reserve(cols);
        inner3.reserve(rows);
        for (int iu = 1; iu < nu; ++iu)
        {
            inner0.push_back(innerAt(iu, 1));
        }
        for (int iv = 1; iv < nv; ++iv)
        {
            inner1.push_back(innerAt(cols, iv));
        }
        for (int iu = nu - 1; iu >= 1; --iu)
        {
            inner2.push_back(innerAt(iu, rows));
        }
        for (int iv = nv - 1; iv >= 1; --iv)
        {
            inner3.push_back(innerAt(1, iv));
        }

        stitchStrip(outer0, inner0);
        stitchStrip(outer1, inner1);
        stitchStrip(outer2, inner2);
        stitchStrip(outer3, inner3);
    }

    std::fclose(f);
    if (outVertexCount)
    {
        *outVertexCount = vertexCount;
    }
    if (outTriangleCount)
    {
        *outTriangleCount = triangleCount;
    }
    return true;
}

int main(int argc, char **argv)
{
    constexpr float kDefaultPixelSpacing = 1.0f;
    constexpr int kDefaultSplitThreshold = 1;
    constexpr int kDefaultSampleSteps = 8;
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
    const pxr::GfVec3f eye =
        camera.found ? camera.worldPosition : (meshCenter + pxr::GfVec3f(0.0f, 0.0f, 5.0f));
    const pxr::GfVec3f lookAt = camera.found ? camera.meshCenter : meshCenter;
    int tmaxComputedEdges = PrecomputePatchEdgeFactors(m,
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
                                                       /*verticalFovDegrees*/ 45.0f);
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
