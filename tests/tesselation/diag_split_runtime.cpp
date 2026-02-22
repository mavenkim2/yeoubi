#include "tesselation/diag_split_runtime.h"
#include "tesselation/diag_split_detail.h"

#include <algorithm>
#include <array>
#include <vector>

namespace ybi::tesselation
{

bool EvaluateLimitPosition(const OpenSubdiv::Far::PatchMap &patchMap,
                           const OpenSubdiv::Far::PatchTable &patchTable,
                           const std::vector<Primvar3> &positions,
                           int ptexFace,
                           float u,
                           float v,
                           pxr::GfVec3f &outPos);

namespace
{

float Lerp(float a, float b, float t)
{
    return a * (1.0f - t) + b * t;
}

void AddTri(TriMesh &out, int i0, int i1, int i2)
{
    out.indices.push_back(i0);
    out.indices.push_back(i1);
    out.indices.push_back(i2);
}

bool AppendEvaluatedPoint(const OpenSubdiv::Far::PatchMap &patchMap,
                          const OpenSubdiv::Far::PatchTable &patchTable,
                          const std::vector<Primvar3> &positions,
                          int ptexFace,
                          float u,
                          float v,
                          TriMesh &out,
                          int &idOut)
{
    pxr::GfVec3f p(0.0f);
    if (!EvaluateLimitPosition(patchMap, patchTable, positions, ptexFace, u, v, p))
    {
        return false;
    }
    idOut = int(out.positions.size());
    out.positions.push_back(p);
    return true;
}

void EdgeUV(const DiagSubPatch &patch, int edge, float t, float &u, float &v)
{
    switch (edge)
    {
        case 0:
            u = Lerp(patch.u0, patch.u1, t);
            v = patch.v0;
            return;
        case 1:
            u = patch.u1;
            v = Lerp(patch.v0, patch.v1, t);
            return;
        case 2:
            u = Lerp(patch.u1, patch.u0, t);
            v = patch.v1;
            return;
        default:
            u = patch.u0;
            v = Lerp(patch.v1, patch.v0, t);
            return;
    }
}

void BuildInnerBoundary(const std::vector<int> &grid, int sideU, int sideV, int edge, std::vector<int> &out)
{
    auto Idx = [sideU, &grid](int u, int v) { return grid[size_t(v * sideU + u)]; };

    out.clear();
    if (edge == 0)
    {
        out.reserve(size_t(sideU));
        for (int u = 0; u < sideU; ++u)
        {
            out.push_back(Idx(u, 0));
        }
        return;
    }
    if (edge == 1)
    {
        out.reserve(size_t(sideV));
        for (int v = 0; v < sideV; ++v)
        {
            out.push_back(Idx(sideU - 1, v));
        }
        return;
    }
    if (edge == 2)
    {
        out.reserve(size_t(sideU));
        for (int u = sideU - 1; u >= 0; --u)
        {
            out.push_back(Idx(u, sideV - 1));
        }
        return;
    }

    out.reserve(size_t(sideV));
    for (int v = sideV - 1; v >= 0; --v)
    {
        out.push_back(Idx(0, v));
    }
}

void StitchEdgeToInner(TriMesh &out, const std::vector<int> &outer, const std::vector<int> &inner)
{
    const int outerSegs = int(outer.size()) - 1;
    const int innerSegs = int(inner.size()) - 1;
    if (outerSegs <= 0 || innerSegs <= 0)
    {
        return;
    }

    int i = 0;
    int j = 0;
    while (i < outerSegs || j < innerSegs)
    {
        const bool advanceOuter =
            (j == innerSegs) ||
            (i < outerSegs && (i == outerSegs - 1 ||
                               float(i + 1) / float(outerSegs) <= float(j + 1) / float(innerSegs)));

        if (advanceOuter)
        {
            AddTri(out, outer[size_t(i)], outer[size_t(i + 1)], inner[size_t(j)]);
            ++i;
        }
        else
        {
            AddTri(out, outer[size_t(i)], inner[size_t(j + 1)], inner[size_t(j)]);
            ++j;
        }
    }
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

        if (faces[pf].fromNgon)
        {
            const int e0 = p.edgeIndex[0];
            if (e0 >= 0 && e0 < int(edges.size()))
            {
                const int r0 = detail::SafeRate(edges, e0, settings.maxRate);
                p.edgeIndex[0] =
                    detail::GetOrCreateChildForNgonBoundary(edges, e0, r0, p.ptexFace, 0.0f, 0.5f);
            }

            const int e3 = p.edgeIndex[3];
            if (e3 >= 0 && e3 < int(edges.size()))
            {
                const int r3 = detail::SafeRate(edges, e3, settings.maxRate);
                p.edgeIndex[3] =
                    detail::GetOrCreateChildForNgonBoundary(edges, e3, r3, p.ptexFace, 1.0f, 0.5f);
            }
        }

        pending.push_back(p);
    }

    while (!pending.empty())
    {
        DiagSubPatch patch = pending.back();
        pending.pop_back();
        out.maxDepthReached = std::max(out.maxDepthReached, patch.depth);

        if (patch.depth >= detail::kDiagSplitMaxDepth || !detail::HasNonUniformEdge(patch, edges))
        {
            out.patches.push_back(patch);
            continue;
        }

        const int r0 = detail::SafeRate(edges, patch.edgeIndex[0], settings.maxRate);
        const int r1 = detail::SafeRate(edges, patch.edgeIndex[1], settings.maxRate);
        const int r2 = detail::SafeRate(edges, patch.edgeIndex[2], settings.maxRate);
        const int r3 = detail::SafeRate(edges, patch.edgeIndex[3], settings.maxRate);
        const bool splitU = (r0 <= 0) || (r2 <= 0);
        const bool splitV = (r1 <= 0) || (r3 <= 0);

        DiagSubPatch a = {};
        DiagSubPatch b = {};

        if (splitU && !splitV)
        {
            detail::SplitPatchU(
                patchMap, patchTable, positions, patch, camera, settings, edges, a, b);
        }
        else if (!splitU && splitV)
        {
            detail::SplitPatchV(
                patchMap, patchTable, positions, patch, camera, settings, edges, a, b);
        }
        else
        {
            const int uWeight = std::max(r0, 1) + std::max(r2, 1);
            const int vWeight = std::max(r1, 1) + std::max(r3, 1);
            if (uWeight >= vWeight)
            {
                detail::SplitPatchU(
                    patchMap, patchTable, positions, patch, camera, settings, edges, a, b);
            }
            else
            {
                detail::SplitPatchV(
                    patchMap, patchTable, positions, patch, camera, settings, edges, a, b);
            }
        }

        pending.push_back(a);
        pending.push_back(b);
        out.splitCount++;
    }

    return out;
}

TriMesh TessellateDiagSplitSubPatches(const OpenSubdiv::Far::PatchMap &patchMap,
                                      const OpenSubdiv::Far::PatchTable &patchTable,
                                      const std::vector<Primvar3> &positions,
                                      const std::vector<DiagSubPatch> &patches,
                                      const std::vector<DiagEdgeState> &edges,
                                      int fallbackEdgeRate,
                                      int &skippedSubPatchesOut)
{
    TriMesh out = {};
    skippedSubPatchesOut = 0;

    std::vector<int> grid;
    std::vector<int> innerBoundary;
    std::array<std::vector<int>, 4> edgeLoops;

    for (const DiagSubPatch &patch : patches)
    {
        const int r0 = detail::PositiveRateOrFallback(edges, patch.edgeIndex[0], fallbackEdgeRate);
        const int r1 = detail::PositiveRateOrFallback(edges, patch.edgeIndex[1], fallbackEdgeRate);
        const int r2 = detail::PositiveRateOrFallback(edges, patch.edgeIndex[2], fallbackEdgeRate);
        const int r3 = detail::PositiveRateOrFallback(edges, patch.edgeIndex[3], fallbackEdgeRate);

        const int Mu = std::max(r0, r2);
        const int Mv = std::max(r1, r3);
        const int sideU = Mu + 1;
        const int sideV = Mv + 1;
        const int base = int(out.positions.size());
        const size_t baseIndices = out.indices.size();

        grid.assign(size_t(sideU * sideV), -1);

        const float insetU = 0.5f / float(std::max(Mu, 2));
        const float insetV = 0.5f / float(std::max(Mv, 2));

        bool ok = true;
        for (int v = 0; v < sideV && ok; ++v)
        {
            for (int u = 0; u < sideU; ++u)
            {
                const float fu = float(u) / float(Mu);
                const float fv = float(v) / float(Mv);
                const float su = Lerp(insetU, 1.0f - insetU, fu);
                const float sv = Lerp(insetV, 1.0f - insetV, fv);
                const float pu = Lerp(patch.u0, patch.u1, su);
                const float pv = Lerp(patch.v0, patch.v1, sv);

                int id = -1;
                if (!AppendEvaluatedPoint(
                        patchMap, patchTable, positions, patch.ptexFace, pu, pv, out, id))
                {
                    ok = false;
                    break;
                }
                grid[size_t(v * sideU + u)] = id;
            }
        }

        std::array<int, 4> corner = {-1, -1, -1, -1};
        if (ok)
        {
            ok &= AppendEvaluatedPoint(
                patchMap, patchTable, positions, patch.ptexFace, patch.u0, patch.v0, out, corner[0]);
            ok &= AppendEvaluatedPoint(
                patchMap, patchTable, positions, patch.ptexFace, patch.u1, patch.v0, out, corner[1]);
            ok &= AppendEvaluatedPoint(
                patchMap, patchTable, positions, patch.ptexFace, patch.u1, patch.v1, out, corner[2]);
            ok &= AppendEvaluatedPoint(
                patchMap, patchTable, positions, patch.ptexFace, patch.u0, patch.v1, out, corner[3]);
        }

        if (ok)
        {
            for (int e = 0; e < 4; ++e)
            {
                const int edgeRate = (e == 0) ? r0 : (e == 1) ? r1 : (e == 2) ? r2 : r3;
                edgeLoops[e].assign(size_t(edgeRate + 1), -1);

                if (e == 0)
                {
                    edgeLoops[e][0] = corner[0];
                    edgeLoops[e][size_t(edgeRate)] = corner[1];
                }
                else if (e == 1)
                {
                    edgeLoops[e][0] = corner[1];
                    edgeLoops[e][size_t(edgeRate)] = corner[2];
                }
                else if (e == 2)
                {
                    edgeLoops[e][0] = corner[2];
                    edgeLoops[e][size_t(edgeRate)] = corner[3];
                }
                else
                {
                    edgeLoops[e][0] = corner[3];
                    edgeLoops[e][size_t(edgeRate)] = corner[0];
                }

                for (int i = 1; i < edgeRate; ++i)
                {
                    const float t = float(i) / float(edgeRate);
                    float pu = 0.0f;
                    float pv = 0.0f;
                    EdgeUV(patch, e, t, pu, pv);
                    int id = -1;
                    if (!AppendEvaluatedPoint(
                            patchMap, patchTable, positions, patch.ptexFace, pu, pv, out, id))
                    {
                        ok = false;
                        break;
                    }
                    edgeLoops[e][size_t(i)] = id;
                }

                if (!ok)
                {
                    break;
                }
            }
        }

        if (!ok)
        {
            out.positions.resize(size_t(base));
            out.indices.resize(baseIndices);
            skippedSubPatchesOut++;
            continue;
        }

        auto Idx = [sideU, &grid](int u, int v) { return grid[size_t(v * sideU + u)]; };
        for (int v = 0; v < Mv; ++v)
        {
            for (int u = 0; u < Mu; ++u)
            {
                const int i0 = Idx(u, v);
                const int i1 = Idx(u + 1, v);
                const int i2 = Idx(u + 1, v + 1);
                const int i3 = Idx(u, v + 1);
                AddTri(out, i0, i1, i2);
                AddTri(out, i0, i2, i3);
            }
        }

        for (int e = 0; e < 4; ++e)
        {
            BuildInnerBoundary(grid, sideU, sideV, e, innerBoundary);
            StitchEdgeToInner(out, edgeLoops[e], innerBoundary);
        }
    }

    return out;
}

} // namespace ybi::tesselation
