#pragma once

#include <vector>

YBI_NAMESPACE_BEGIN

struct float3;
struct Device;

enum BVHFlags : int
{
    USE_CLUSTERS = (1u << 0u),
};

struct BVH
{
    BVHFlags flags;

    BVH();
};

struct Mesh
{
    const float3 *positions;
    const int *indices;
    uint32_t numVertices;
    uint32_t numIndices;

    Mesh(float3 *positions, int *indices, uint32_t numVertices, uint32_t numIndices);
};

struct Scene
{
    BVH bvh;
    std::vector<Mesh> meshes;
    Device *device;
};

YBI_NAMESPACE_END
