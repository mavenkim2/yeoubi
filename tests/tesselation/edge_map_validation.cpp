#include "tesselation/edge_map_validation.h"

#include <algorithm>
#include <unordered_map>

EdgeMapChecks RunEdgeMapChecks(const SelectedSubdivMesh &m,
                               const std::vector<SubdivisionPatch> &patches,
                               const SubdivisionEdgeMap &edgeMap)
{
    EdgeMapChecks checks = {};

    // 1) Every coarse edge touched by a non-quad face must have a midpoint vertex allocated.
    std::unordered_map<int, size_t> faceStart;
    faceStart.reserve(m.faceVertexCounts.size());
    size_t cursor = 0;
    for (size_t f = 0; f < m.faceVertexCounts.size(); ++f)
    {
        faceStart[int(f)] = cursor;
        cursor += std::max(0, m.faceVertexCounts[f]);
    }

    std::unordered_map<uint64_t, int> checkedNgonEdgeKeys;
    checkedNgonEdgeKeys.reserve(m.faceVertexIndices.size());
    for (size_t f = 0; f < m.faceVertexCounts.size(); ++f)
    {
        const int n = m.faceVertexCounts[f];
        if (n <= 4)
        {
            continue;
        }
        const size_t start = faceStart[int(f)];
        if (start + size_t(n) > m.faceVertexIndices.size())
        {
            continue;
        }
        for (int i = 0; i < n; ++i)
        {
            const int v0 = m.faceVertexIndices[start + i];
            const int v1 = m.faceVertexIndices[start + ((i + 1) % n)];
            const uint64_t key = MakeEdgeKey(v0, v1);
            if (!checkedNgonEdgeKeys.insert({key, 1}).second)
            {
                continue;
            }
            const auto it = edgeMap.find(key);
            if (it == edgeMap.end() || it->second.midpointVertex < 0)
            {
                checks.missingNgonMidpoints++;
            }
        }
    }

    // 2) A midpoint vertex id must map to exactly one coarse edge.
    std::unordered_map<int, uint64_t> midpointToEdge;
    midpointToEdge.reserve(edgeMap.size());
    for (const auto &kv : edgeMap)
    {
        const uint64_t edgeKey = kv.first;
        const SubdivisionEdge &edge = kv.second;
        if (edge.midpointVertex < 0)
        {
            continue;
        }
        auto it = midpointToEdge.find(edge.midpointVertex);
        if (it == midpointToEdge.end())
        {
            midpointToEdge.insert({edge.midpointVertex, edgeKey});
        }
        else if (it->second != edgeKey)
        {
            checks.duplicateMidpointVertices++;
        }
    }

    // 3) Generated edges used by non-quad patches must exist in edgeMap.
    for (const SubdivisionPatch &patch : patches)
    {
        if (patch.coarseFace < 0 || patch.coarseFace >= int(m.faceVertexCounts.size()))
        {
            continue;
        }
        if (m.faceVertexCounts[patch.coarseFace] <= 4)
        {
            continue;
        }
        const uint64_t edgeKeys[4] = {MakeEdgeKey(patch.verts[0], patch.verts[1]),
                                      MakeEdgeKey(patch.verts[1], patch.verts[2]),
                                      MakeEdgeKey(patch.verts[2], patch.verts[3]),
                                      MakeEdgeKey(patch.verts[3], patch.verts[0])};
        for (int i = 0; i < 4; ++i)
        {
            if (edgeMap.find(edgeKeys[i]) == edgeMap.end())
            {
                checks.missingGeneratedPatchEdges++;
            }
        }
    }

    // 4) Stored flags must match faceCount-derived flags.
    for (const auto &kv : edgeMap)
    {
        const SubdivisionEdge &edge = kv.second;
        const bool expectedBoundary = (edge.faceCount == 1);
        const bool expectedNonManifold = (edge.faceCount > 2);
        if (edge.boundary != expectedBoundary)
        {
            checks.badBoundaryFlags++;
        }
        if (edge.nonManifold != expectedNonManifold)
        {
            checks.badNonManifoldFlags++;
        }
    }

    checks.ok = (checks.missingNgonMidpoints == 0) && (checks.duplicateMidpointVertices == 0) &&
                (checks.missingGeneratedPatchEdges == 0) && (checks.badBoundaryFlags == 0) &&
                (checks.badNonManifoldFlags == 0);
    return checks;
}
