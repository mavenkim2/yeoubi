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
};

struct Mesh
{
    float3 *positions;
    int *indices;

    Mesh(float3 *positions, int *indices);
};

struct Scene
{
    BVH bvh;
    std::vector<Mesh> meshes;
    Device *device;
};

YBI_NAMESPACE_END
