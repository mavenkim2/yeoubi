#include "tesselation/diag_split_runtime.h"

#include "util/assert.h"

#include <algorithm>
#include <utility>

namespace ybi::tesselation
{

bool EvaluateLimitPosition(const OpenSubdiv::Far::PatchMap &patchMap,
                           const OpenSubdiv::Far::PatchTable &patchTable,
                           const std::vector<Primvar3> &positions,
                           int ptexFace,
                           float u,
                           float v,
                           pxr::GfVec3f &outPos);

int ComputeDiagSplitEdgeFactor(const OpenSubdiv::Far::PatchMap &patchMap,
                               const OpenSubdiv::Far::PatchTable &patchTable,
                               const std::vector<Primvar3> &positions,
                               int ptexFace,
                               const pxr::GfVec2f &uvStart,
                               const pxr::GfVec2f &uvEnd,
                               const EdgeFactorCamera &camera,
                               const EdgeFactorSettings &settings);

namespace
{

constexpr int kDiagSplitMaxDepth = 16;

bool HasNonUniformEdge(const DiagSubPatch &patch, const std::vector<DiagEdgeState> &edges)
{
    for (int e = 0; e < 4; ++e)
    {
        const int edgeId = patch.edgeIndex[e];
        if (edgeId < 0 || edgeId >= int(edges.size()) || edges[edgeId].rate <= 0)
        {
            return true;
        }
    }
    return false;
}

int SafeRate(const std::vector<DiagEdgeState> &edges, int edgeId, int fallback)
{
    if (edgeId < 0 || edgeId >= int(edges.size()))
    {
        return fallback;
    }
    return edges[edgeId].rate;
}

int PositiveRateOrFallback(const std::vector<DiagEdgeState> &edges, int edgeId, int fallback)
{
    const int r = SafeRate(edges, edgeId, fallback);
    if (r <= 0)
    {
        return std::max(1, fallback);
    }
    return r;
}

int AppendEdgeRate(std::vector<DiagEdgeState> &edges, int rate)
{
    edges.push_back(DiagEdgeState{-1, -1, rate});
    return int(edges.size()) - 1;
}

std::pair<int, int> GetOrCreateKnownChildren(std::vector<DiagEdgeState> &edges,
                                             int parentEdgeId,
                                             int parentRate)
{
    DiagEdgeState &parent = edges[parentEdgeId];
    if (parent.child0 >= 0 && parent.child1 >= 0)
    {
        return {parent.child0, parent.child1};
    }
    const int k = std::max(1, parentRate / 2);
    parent.child0 = AppendEdgeRate(edges, k);
    parent.child1 = AppendEdgeRate(edges, std::max(1, parentRate - k));
    return {parent.child0, parent.child1};
}

std::pair<int, int> GetOrCreateComputedChildren(const OpenSubdiv::Far::PatchMap &patchMap,
                                                const OpenSubdiv::Far::PatchTable &patchTable,
                                                const std::vector<Primvar3> &positions,
                                                int ptexFace,
                                                int parentEdgeId,
                                                const pxr::GfVec2f &uvStart,
                                                const pxr::GfVec2f &uvMid,
                                                const pxr::GfVec2f &uvEnd,
                                                const EdgeFactorCamera &camera,
                                                const EdgeFactorSettings &settings,
                                                std::vector<DiagEdgeState> &edges)
{
    DiagEdgeState &parent = edges[parentEdgeId];
    if (parent.child0 >= 0 && parent.child1 >= 0)
    {
        return {parent.child0, parent.child1};
    }
    const int r0 = ComputeDiagSplitEdgeFactor(
        patchMap, patchTable, positions, ptexFace, uvStart, uvMid, camera, settings);
    const int r1 = ComputeDiagSplitEdgeFactor(
        patchMap, patchTable, positions, ptexFace, uvMid, uvEnd, camera, settings);
    parent.child0 = AppendEdgeRate(edges, r0);
    parent.child1 = AppendEdgeRate(edges, r1);
    return {parent.child0, parent.child1};
}

void SplitPatchU(const OpenSubdiv::Far::PatchMap &patchMap,
                 const OpenSubdiv::Far::PatchTable &patchTable,
                 const std::vector<Primvar3> &positions,
                 const DiagSubPatch &patch,
                 const EdgeFactorCamera &camera,
                 const EdgeFactorSettings &settings,
                 std::vector<DiagEdgeState> &edges,
                 DiagSubPatch &a,
                 DiagSubPatch &b)
{
    const int r0 = SafeRate(edges, patch.edgeIndex[0], settings.maxRate);
    const int r2 = SafeRate(edges, patch.edgeIndex[2], settings.maxRate);
    const bool known0 = (r0 > 1);
    const bool known2 = (r2 > 1);

    float splitU01 = 0.5f;
    if (known0 && !known2)
    {
        splitU01 = float(r0 / 2) / float(r0);
    }
    else if (!known0 && known2)
    {
        splitU01 = 1.0f - (float(r2 / 2) / float(r2));
    }
    const float um = patch.u0 * (1.0f - splitU01) + patch.u1 * splitU01;

    a = patch;
    b = patch;
    a.u1 = um;
    b.u0 = um;
    a.depth = patch.depth + 1;
    b.depth = patch.depth + 1;

    const pxr::GfVec2f uvBottom0(patch.u0, patch.v0);
    const pxr::GfVec2f uvBottomM(um, patch.v0);
    const pxr::GfVec2f uvBottom1(patch.u1, patch.v0);
    const pxr::GfVec2f uvTop0(patch.u0, patch.v1);
    const pxr::GfVec2f uvTopM(um, patch.v1);
    const pxr::GfVec2f uvTop1(patch.u1, patch.v1);
    const pxr::GfVec2f uvCenter0(um, patch.v0);
    const pxr::GfVec2f uvCenter1(um, patch.v1);

    YBI_ASSERT(patch.edgeIndex[0] >= 0 && patch.edgeIndex[0] < int(edges.size()));
    YBI_ASSERT(patch.edgeIndex[2] >= 0 && patch.edgeIndex[2] < int(edges.size()));

    if (known0)
    {
        const auto children = GetOrCreateKnownChildren(edges, patch.edgeIndex[0], r0);
        a.edgeIndex[0] = children.first;
        b.edgeIndex[0] = children.second;
    }
    else
    {
        const auto children = GetOrCreateComputedChildren(patchMap,
                                                          patchTable,
                                                          positions,
                                                          patch.ptexFace,
                                                          patch.edgeIndex[0],
                                                          uvBottom0,
                                                          uvBottomM,
                                                          uvBottom1,
                                                          camera,
                                                          settings,
                                                          edges);
        a.edgeIndex[0] = children.first;
        b.edgeIndex[0] = children.second;
    }

    if (known2)
    {
        const auto children = GetOrCreateKnownChildren(edges, patch.edgeIndex[2], r2);
        b.edgeIndex[2] = children.first;
        a.edgeIndex[2] = children.second;
    }
    else
    {
        const auto children = GetOrCreateComputedChildren(patchMap,
                                                          patchTable,
                                                          positions,
                                                          patch.ptexFace,
                                                          patch.edgeIndex[2],
                                                          uvTop1,
                                                          uvTopM,
                                                          uvTop0,
                                                          camera,
                                                          settings,
                                                          edges);
        b.edgeIndex[2] = children.first;
        a.edgeIndex[2] = children.second;
    }

    a.edgeIndex[3] = patch.edgeIndex[3];
    b.edgeIndex[1] = patch.edgeIndex[1];
    const int centerRate = ComputeDiagSplitEdgeFactor(
        patchMap, patchTable, positions, patch.ptexFace, uvCenter0, uvCenter1, camera, settings);
    const int centerEdge = AppendEdgeRate(edges, centerRate);
    a.edgeIndex[1] = centerEdge;
    b.edgeIndex[3] = centerEdge;
}

void SplitPatchV(const OpenSubdiv::Far::PatchMap &patchMap,
                 const OpenSubdiv::Far::PatchTable &patchTable,
                 const std::vector<Primvar3> &positions,
                 const DiagSubPatch &patch,
                 const EdgeFactorCamera &camera,
                 const EdgeFactorSettings &settings,
                 std::vector<DiagEdgeState> &edges,
                 DiagSubPatch &a,
                 DiagSubPatch &b)
{
    const int r1 = SafeRate(edges, patch.edgeIndex[1], settings.maxRate);
    const int r3 = SafeRate(edges, patch.edgeIndex[3], settings.maxRate);
    const bool known1 = (r1 > 1);
    const bool known3 = (r3 > 1);

    float splitV01 = 0.5f;
    if (known1 && !known3)
    {
        splitV01 = float(r1 / 2) / float(r1);
    }
    else if (!known1 && known3)
    {
        splitV01 = 1.0f - (float(r3 / 2) / float(r3));
    }
    const float vm = patch.v0 * (1.0f - splitV01) + patch.v1 * splitV01;

    a = patch;
    b = patch;
    a.v1 = vm;
    b.v0 = vm;
    a.depth = patch.depth + 1;
    b.depth = patch.depth + 1;

    const pxr::GfVec2f uvRight0(patch.u1, patch.v0);
    const pxr::GfVec2f uvRightM(patch.u1, vm);
    const pxr::GfVec2f uvRight1(patch.u1, patch.v1);
    const pxr::GfVec2f uvLeft0(patch.u0, patch.v0);
    const pxr::GfVec2f uvLeftM(patch.u0, vm);
    const pxr::GfVec2f uvLeft1(patch.u0, patch.v1);
    const pxr::GfVec2f uvCenter0(patch.u0, vm);
    const pxr::GfVec2f uvCenter1(patch.u1, vm);

    YBI_ASSERT(patch.edgeIndex[1] >= 0 && patch.edgeIndex[1] < int(edges.size()));
    YBI_ASSERT(patch.edgeIndex[3] >= 0 && patch.edgeIndex[3] < int(edges.size()));

    if (known1)
    {
        const auto children = GetOrCreateKnownChildren(edges, patch.edgeIndex[1], r1);
        a.edgeIndex[1] = children.first;
        b.edgeIndex[1] = children.second;
    }
    else
    {
        const auto children = GetOrCreateComputedChildren(patchMap,
                                                          patchTable,
                                                          positions,
                                                          patch.ptexFace,
                                                          patch.edgeIndex[1],
                                                          uvRight0,
                                                          uvRightM,
                                                          uvRight1,
                                                          camera,
                                                          settings,
                                                          edges);
        a.edgeIndex[1] = children.first;
        b.edgeIndex[1] = children.second;
    }

    if (known3)
    {
        const auto children = GetOrCreateKnownChildren(edges, patch.edgeIndex[3], r3);
        b.edgeIndex[3] = children.first;
        a.edgeIndex[3] = children.second;
    }
    else
    {
        const auto children = GetOrCreateComputedChildren(patchMap,
                                                          patchTable,
                                                          positions,
                                                          patch.ptexFace,
                                                          patch.edgeIndex[3],
                                                          uvLeft1,
                                                          uvLeftM,
                                                          uvLeft0,
                                                          camera,
                                                          settings,
                                                          edges);
        b.edgeIndex[3] = children.first;
        a.edgeIndex[3] = children.second;
    }

    a.edgeIndex[0] = patch.edgeIndex[0];
    b.edgeIndex[2] = patch.edgeIndex[2];
    const int centerRate = ComputeDiagSplitEdgeFactor(
        patchMap, patchTable, positions, patch.ptexFace, uvCenter0, uvCenter1, camera, settings);
    const int centerEdge = AppendEdgeRate(edges, centerRate);
    a.edgeIndex[2] = centerEdge;
    b.edgeIndex[0] = centerEdge;
}

} // namespace

DiagSplitBuildResult BuildDiagSplitSubPatches(const OpenSubdiv::Far::PatchMap &patchMap,
                                              const OpenSubdiv::Far::PatchTable &patchTable,
                                              const std::vector<Primvar3> &positions,
                                              const std::vector<PtexFaceAdj> &faces,
                                              int numPtexFaces,
                                              const EdgeFactorCamera &camera,
                                              const EdgeFactorSettings &settings,
                                              std::vector<DiagEdgeState> &edges)
{
    DiagSplitBuildResult out = {};
    std::vector<DiagSubPatch> pending;
    pending.reserve(size_t(numPtexFaces));
    for (int pf = 0; pf < numPtexFaces; ++pf)
    {
        if (!patchMap.FindPatch(pf, 0.5f, 0.5f))
        {
            out.skippedPtexFaces++;
            continue;
        }
        DiagSubPatch p = {};
        p.ptexFace = pf;
        p.edgeIndex = faces[pf].edgeIndex;
        pending.push_back(p);
    }

    while (!pending.empty())
    {
        DiagSubPatch patch = pending.back();
        pending.pop_back();
        out.maxDepthReached = std::max(out.maxDepthReached, patch.depth);
        if (patch.depth >= kDiagSplitMaxDepth || !HasNonUniformEdge(patch, edges))
        {
            out.patches.push_back(patch);
            continue;
        }

        const int r0 = SafeRate(edges, patch.edgeIndex[0], settings.maxRate);
        const int r1 = SafeRate(edges, patch.edgeIndex[1], settings.maxRate);
        const int r2 = SafeRate(edges, patch.edgeIndex[2], settings.maxRate);
        const int r3 = SafeRate(edges, patch.edgeIndex[3], settings.maxRate);
        const bool splitU = (r0 <= 0) || (r2 <= 0);
        const bool splitV = (r1 <= 0) || (r3 <= 0);

        DiagSubPatch a = {};
        DiagSubPatch b = {};
        if (splitU && !splitV)
        {
            SplitPatchU(patchMap, patchTable, positions, patch, camera, settings, edges, a, b);
        }
        else if (!splitU && splitV)
        {
            SplitPatchV(patchMap, patchTable, positions, patch, camera, settings, edges, a, b);
        }
        else
        {
            const int uWeight = std::max(r0, 1) + std::max(r2, 1);
            const int vWeight = std::max(r1, 1) + std::max(r3, 1);
            if (uWeight >= vWeight)
            {
                SplitPatchU(
                    patchMap, patchTable, positions, patch, camera, settings, edges, a, b);
            }
            else
            {
                SplitPatchV(
                    patchMap, patchTable, positions, patch, camera, settings, edges, a, b);
            }
        }

        pending.push_back(a);
        pending.push_back(b);
        out.splitCount++;
    }
    return out;
}

std::vector<DiagEdgeState> BuildDiagEdgesFromRates(const std::vector<int> &edgeRates)
{
    std::vector<DiagEdgeState> edges;
    edges.reserve(edgeRates.size());
    for (int rate : edgeRates)
    {
        edges.push_back(DiagEdgeState{-1, -1, rate});
    }
    return edges;
}

TriMesh TessellateDiagSplitNoStitch(const OpenSubdiv::Far::PatchMap &patchMap,
                                    const OpenSubdiv::Far::PatchTable &patchTable,
                                    const std::vector<Primvar3> &positions,
                                    const std::vector<DiagSubPatch> &patches,
                                    const std::vector<DiagEdgeState> &edges,
                                    int fallbackEdgeRate,
                                    int &skippedSubPatchesOut)
{
    TriMesh out = {};
    skippedSubPatchesOut = 0;
    for (const DiagSubPatch &patch : patches)
    {
        const int r0 = PositiveRateOrFallback(edges, patch.edgeIndex[0], fallbackEdgeRate);
        const int r1 = PositiveRateOrFallback(edges, patch.edgeIndex[1], fallbackEdgeRate);
        const int r2 = PositiveRateOrFallback(edges, patch.edgeIndex[2], fallbackEdgeRate);
        const int r3 = PositiveRateOrFallback(edges, patch.edgeIndex[3], fallbackEdgeRate);
        const int Mu = std::max(r0, r2);
        const int Mv = std::max(r1, r3);
        const int sideU = Mu + 1;
        const int sideV = Mv + 1;
        const int base = int(out.positions.size());

        bool ok = true;
        for (int v = 0; v < sideV && ok; ++v)
        {
            for (int u = 0; u < sideU; ++u)
            {
                const float fu = float(u) / float(Mu);
                const float fv = float(v) / float(Mv);
                const float pu = patch.u0 * (1.0f - fu) + patch.u1 * fu;
                const float pv = patch.v0 * (1.0f - fv) + patch.v1 * fv;
                pxr::GfVec3f p(0.0f);
                if (!EvaluateLimitPosition(patchMap, patchTable, positions, patch.ptexFace, pu, pv, p))
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
            skippedSubPatchesOut++;
            continue;
        }

        auto Idx = [sideU, base](int u, int v) { return base + v * sideU + u; };
        for (int v = 0; v < Mv; ++v)
        {
            for (int u = 0; u < Mu; ++u)
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

} // namespace ybi::tesselation
