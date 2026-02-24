#pragma once

#include <algorithm>
#include <cstdint>
#include <unordered_map>

struct SubdivisionEdge
{
    int v0 = -1;
    int v1 = -1;
    int midpointVertex = -1;
    int faceCount = 0;
    int firstFace = -1;
    int secondFace = -1;
    bool boundary = false;
    bool nonManifold = false;
};

using SubdivisionEdgeMap = std::unordered_map<uint64_t, SubdivisionEdge>;

struct SubdivisionPatch
{
    int verts[4] = {-1, -1, -1, -1};
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
