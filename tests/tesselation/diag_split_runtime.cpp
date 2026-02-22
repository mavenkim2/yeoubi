#include "tesselation/diag_split_runtime.h"
#include "tesselation/diag_split_detail.h"

#include <algorithm>

namespace ybi::tesselation
{

bool EvaluateLimitPosition(const OpenSubdiv::Far::PatchMap &patchMap,
                           const OpenSubdiv::Far::PatchTable &patchTable,
                           const std::vector<Primvar3> &positions,
                           int ptexFace,
                           float u,
                           float v,
                           pxr::GfVec3f &outPos);

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

    for (const DiagSubPatch &patch : patches)
    {
        const int r0 = detail::PositiveRateOrFallback(edges, patch.edgeIndex[0], fallbackEdgeRate);
        const int r1 = detail::PositiveRateOrFallback(edges, patch.edgeIndex[1], fallbackEdgeRate);
        const int r2 = detail::PositiveRateOrFallback(edges, patch.edgeIndex[2], fallbackEdgeRate);
        const int r3 = detail::PositiveRateOrFallback(edges, patch.edgeIndex[3], fallbackEdgeRate);

        const int Mu = std::max(r0, r2);
        const int Mv = std::max(r1, r3);
        if (Mu < 2 || Mv < 2)
        {
            skippedSubPatchesOut++;
            continue;
        }

        const int sideU = Mu + 1;
        const int sideV = Mv + 1;
        const int base = int(out.positions.size());
        const float insetU = 0.5f / float(Mu);
        const float insetV = 0.5f / float(Mv);

        bool ok = true;
        for (int v = 0; v < sideV && ok; ++v)
        {
            for (int u = 0; u < sideU; ++u)
            {
                const float fu = float(u) / float(Mu);
                const float fv = float(v) / float(Mv);
                const float su = insetU + fu * (1.0f - 2.0f * insetU);
                const float sv = insetV + fv * (1.0f - 2.0f * insetV);
                const float pu = patch.u0 * (1.0f - su) + patch.u1 * su;
                const float pv = patch.v0 * (1.0f - sv) + patch.v1 * sv;

                pxr::GfVec3f p(0.0f);
                if (!EvaluateLimitPosition(
                        patchMap, patchTable, positions, patch.ptexFace, pu, pv, p))
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
