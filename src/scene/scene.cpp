#include "scene.h"
#include "util/float3.h"

YBI_NAMESPACE_BEGIN

BVH::BVH() : flags(BVHFlags(0)) {}

Mesh::Mesh(Array<float3> &&pos, Array<int> &&idx)
    : positions(std::move(pos)), indices(std::move(idx))
{
}

Curves::Curves(float3 *positions, int *curveVertexOffsets, int numVertices, int numCurves)
    : positions(positions), curveVertexOffsets(curveVertexOffsets), numVertices(numVertices),
      numCurves(numCurves)
{
}

YBI_NAMESPACE_END
