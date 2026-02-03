#include "scene.h"
#include "../util/float3.h"

YBI_NAMESPACE_BEGIN

BVH::BVH() : flags(BVHFlags(0)) {}
Mesh::Mesh(float3 *positions, int *indices, uint32_t numVertices, uint32_t numIndices)
    : positions(positions), indices(indices), numVertices(numVertices), numIndices(numIndices)
{
}

YBI_NAMESPACE_END
