#pragma once

#include <vector>

YBI_NAMESPACE_BEGIN

struct float3;
struct Device;

enum BVHFlags : int
{
    USE_CLUSTERS = (1u << 0u),
};

enum CurveFlags : int
{
    CURVE_FLAGS_RIBBON = (1u << 0u),
    CURVE_FLAGS_TUBE   = (1u << 1u),

    CURVE_FLAGS_LINEAR = (1u << 2u),
    CURVE_FLAGS_CUBIC  = (1u << 3u),

    CURVE_FLAGS_BSPLINE = (1u << 4u),
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

struct Curves
{
    // num_offsets = num_curves + 1
    int *curveVertexOffsets;
    int numCurves;

    int curveFlags;
};

struct Scene
{
    BVH bvh;
    std::vector<Mesh> meshes;
    Device *device;
};

YBI_NAMESPACE_END
