#pragma once

#include <algorithm>
#include <cstdint>
#include <pxr/base/gf/vec2f.h>
#include <unordered_map>

inline constexpr int SUBDIV_EDGE_FACTOR_NON_UNIFORM = -1;
inline constexpr int SUBDIV_EDGE_FACTOR_UNINITIALIZED = -2;

struct SubdivisionEdge
{
    int v0 = -1;
    int v1 = -1;
    int sampleVStart = -1;
    int sampleVEnd = -1;
    int midpointVertex = -1;
    int edgeVertexIndexStart = -1;
    int storedPtexFaceId = -1;
    pxr::GfVec2f storedUv0 = pxr::GfVec2f(0.0f, 0.0f);
    pxr::GfVec2f storedUv1 = pxr::GfVec2f(0.0f, 0.0f);
    bool hasStoredPatchParams = false;
    bool transitionedUninitializedToUniform = false;
    int tmaxEdgeFactor = SUBDIV_EDGE_FACTOR_UNINITIALIZED;
    int faceCount = 0;
    int firstFace = -1;
    int secondFace = -1;
    bool boundary = false;
    bool nonManifold = false;

    bool split = false;
};

using SubdivisionEdgeMap = std::unordered_map<uint64_t, SubdivisionEdge>;

struct SubdivisionPatch
{
    int verts[4] = {-1, -1, -1, -1};
    pxr::GfVec2f uv[4];
    int coarseFace = -1;
    int quadrant = 0;
    int ptexFaceId = -1;
};

inline uint64_t MakeEdgeKey(int v0, int v1)
{
    const int lo = std::min(v0, v1);
    const int hi = std::max(v0, v1);
    return (uint64_t(uint32_t(lo)) << 32) | uint64_t(uint32_t(hi));
}
