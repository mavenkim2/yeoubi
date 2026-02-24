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

static SubdivisionEdge &GetOrCreateEdge(SubdivisionEdgeMap &edgeMap, int v0, int v1)
{
    const uint64_t key = MakeEdgeKey(v0, v1);
    auto it = edgeMap.find(key);
    if (it == edgeMap.end())
    {
        SubdivisionEdge edge = {};
        edge.v0 = v0;
        edge.v1 = v1;
        it = edgeMap.emplace(key, edge).first;
    }
    return it->second;
}

static int ComputeUniformSplitPosStored(const SubdivisionEdge &edge)
{
    YBI_ASSERT(edge.tmaxEdgeFactor >= 2);
    YBI_ASSERT(edge.edgeVertexIndexStart >= 0);
    YBI_ASSERT(edge.midpointVertex >= 0);
    const int splitPos = edge.midpointVertex - edge.edgeVertexIndexStart + 1;
    YBI_ASSERT(splitPos >= 1 && splitPos < edge.tmaxEdgeFactor);
    return splitPos;
}

static void InitializeUniformChildEdgeRangeFromParent(const SubdivisionEdge &parent,
                                                      int splitPosStored,
                                                      int childFactor,
                                                      SubdivisionEdge *child)
{
    YBI_ASSERT(child);
    YBI_ASSERT(parent.tmaxEdgeFactor >= 2);
    YBI_ASSERT(parent.edgeVertexIndexStart >= 0);

    const int p0 = parent.v0;
    const int p1 = parent.v1;
    const int mid = parent.midpointVertex;

    const bool usesLeft =
        (child->v0 == p0 && child->v1 == mid) || (child->v0 == mid && child->v1 == p0);
    const bool usesRight =
        (child->v0 == mid && child->v1 == p1) || (child->v0 == p1 && child->v1 == mid);
    YBI_ASSERT(usesLeft || usesRight);
    YBI_ASSERT(childFactor >= 1);

    const int intervalInternalCount = childFactor - 1;
    int expectedStart = -1;
    if (intervalInternalCount > 0)
    {
        const int intervalBase =
            usesLeft ? parent.edgeVertexIndexStart : (parent.edgeVertexIndexStart + splitPosStored);
        const bool forward =
            usesLeft ? (child->v0 == p0 && child->v1 == mid) : (child->v0 == mid && child->v1 == p1);
        expectedStart = forward ? intervalBase : (intervalBase + intervalInternalCount - 1);
    }

    if (child->edgeVertexIndexStart < 0)
    {
        child->edgeVertexIndexStart = expectedStart;
    }
    else if (child->edgeVertexIndexStart != expectedStart && expectedStart >= 0)
    {
        // Keep existing range initialization if another patch already established it.
    }

    if (childFactor >= 2)
    {
        if (child->edgeVertexIndexStart < 0)
        {
            child->edgeVertexIndexStart = expectedStart;
        }
        YBI_ASSERT(child->edgeVertexIndexStart >= 0);
        child->midpointVertex = child->edgeVertexIndexStart + (childFactor / 2 - 1);
    }
    else
    {
        child->midpointVertex = -1;
    }
}

static int
EnsurePatchEdgeFactor(const SelectedSubdivMesh &m,
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
    SubdivisionEdge &edge = GetOrCreateEdge(edgeMap, patch.verts[edgeIndex], patch.verts[next]);
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
        const bool coarseFaceIsQuad = (patch.coarseFace >= 0) &&
                                      (patch.coarseFace < int(m.faceVertexCounts.size())) &&
                                      (m.faceVertexCounts[patch.coarseFace] == 4);
        if (allowQuadCoarseUpgrade && coarseFaceIsQuad && edge.tmaxEdgeFactor == 1)
        {
            edge.tmaxEdgeFactor = 2;
        }
        if (edge.tmaxEdgeFactor >= 1)
        {
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
    return edge.tmaxEdgeFactor;
}

static int PrecomputePatchEdgeFactors(
    const SelectedSubdivMesh &m,
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
            SubdivisionEdge &opp = GetOrCreateEdge(edgeMap, ov0, ov1);
            opp.tmaxEdgeFactor = 2;
            if (opp.edgeVertexIndexStart < 0)
            {
                opp.edgeVertexIndexStart = nextGeneratedVertexId++;
            }
            opp.midpointVertex = opp.edgeVertexIndexStart;
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
            SubdivisionEdge &edge = GetOrCreateEdge(edgeMap, v0, v1);

            int mid = edge.midpointVertex;
            int t = edgeFactor[e];
            float alpha = 0.5f;
            int leftFactor = SUBDIV_EDGE_FACTOR_UNINITIALIZED;
            int rightFactor = SUBDIV_EDGE_FACTOR_UNINITIALIZED;
            if (t != SUBDIV_EDGE_FACTOR_NON_UNIFORM)
            {
                YBI_ASSERT(t > 1);
                if (edge.edgeVertexIndexStart < 0 || mid < 0)
                {
                    const int splitPosStored = t / 2;
                    edge.edgeVertexIndexStart = nextGeneratedVertexId;
                    nextGeneratedVertexId += (t - 1);
                    edge.midpointVertex = edge.edgeVertexIndexStart + (splitPosStored - 1);
                    mid = edge.midpointVertex;
                }
                const int splitPosStored = ComputeUniformSplitPosStored(edge);
                const bool localForward = (edge.v0 == v0 && edge.v1 == v1);
                leftFactor = localForward ? splitPosStored : (t - splitPosStored);
                rightFactor = t - leftFactor;
                YBI_ASSERT(leftFactor + rightFactor == t);
                alpha = float(leftFactor) / float(t);
            }

            if (mid < 0)
            {
                mid = nextGeneratedVertexId++;
                edge.midpointVertex = mid;
            }
            splitVerts[i] = mid;
            splitUV[i] = patch.uv[e] * (1.0f - alpha) + patch.uv[n] * alpha;

            SubdivisionEdge &edgeA = GetOrCreateEdge(edgeMap, v0, mid);
            SubdivisionEdge &edgeB = GetOrCreateEdge(edgeMap, mid, v1);
            edgeA.tmaxEdgeFactor = leftFactor;
            edgeB.tmaxEdgeFactor = rightFactor;
            edgeA.boundary = edge.boundary;
            edgeB.boundary = edge.boundary;
            if (t != SUBDIV_EDGE_FACTOR_NON_UNIFORM)
            {
                const int splitPosStored = ComputeUniformSplitPosStored(edge);
                InitializeUniformChildEdgeRangeFromParent(edge, splitPosStored, leftFactor, &edgeA);
                InitializeUniformChildEdgeRangeFromParent(
                    edge, splitPosStored, rightFactor, &edgeB);
            }
            else
            {
                edgeA.edgeVertexIndexStart = -1;
                edgeB.edgeVertexIndexStart = -1;
            }
        }

        SubdivisionEdge &seam = GetOrCreateEdge(edgeMap, splitVerts[0], splitVerts[1]);
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

static int
GetOrAllocateMidpointVertex(SubdivisionEdgeMap &edgeMap, int v0, int v1, int &nextVertexId)
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
                const int midNext = GetOrAllocateMidpointVertex(edgeMap, v, vNext, nextVertexId);
                const int midPrev = GetOrAllocateMidpointVertex(edgeMap, vPrev, v, nextVertexId);

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

int main(int argc, char **argv)
{
    if (argc < 2)
    {
        std::fprintf(stderr, "Usage: %s <selected-subdiv.json> [level>=1]\n", argv[0]);
        return 2;
    }

    const std::string inJson = argv[1];
    int level = 1;
    if (argc >= 3)
    {
        level = std::max(1, std::atoi(argv[2]));
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
                                                       /*sampleStepsN*/ 8,
                                                       /*targetPixelSpacing*/ 0.25f,
                                                       /*splitThreshold*/ 2,
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
                         /*sampleStepsN*/ 8,
                         /*targetPixelSpacing*/ 0.25f,
                         /*splitThreshold*/ 2,
                         eye,
                         lookAt,
                         /*viewportWidth*/ 1920,
                         /*viewportHeight*/ 1080,
                         /*verticalFovDegrees*/ 45.0f,
                         &tmaxComputedEdges);

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
    std::printf("  edgeMapChecks=ok\n");
    std::printf("  controlCageUniqueEdges=%zu boundaryEdges=%d\n",
                edgeMap.size(),
                CountBoundaryEdges(edgeMap));
    std::printf("  controlCageEdgesWithOver2Faces=%d\n", edgesWithOver2Faces);

    delete patchTable;
    delete refiner;
    return 0;
}
