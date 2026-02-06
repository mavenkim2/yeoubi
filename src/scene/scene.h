#pragma once

#include "util/array.h"
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
    CURVE_FLAGS_TUBE = (1u << 1u),

    CURVE_FLAGS_LINEAR = (1u << 2u),
    CURVE_FLAGS_CUBIC = (1u << 3u),

    CURVE_FLAGS_BSPLINE = (1u << 4u),
};

struct BVH
{
    BVHFlags flags;

    BVH();
};

struct Mesh
{
    Array<float3> positions;
    Array<int> indices;

    Mesh() = default;
    ~Mesh() = default;
    Mesh(Mesh &&other) = default;

    Mesh(Array<float3> &&pos, Array<int> &&idx);
};

struct Curves
{
    const float3 *positions;

    // num_offsets = num_curves + 1
    int *curveVertexOffsets;
    int numCurves;
    int numVertices;

    int curveFlags;

    Curves() = default;
    ~Curves() = default;

    Curves(float3 *positions, int *curveVertexOffsets, int numVertices, int numCurves);
};

struct Scene
{
    BVH bvh;
    // Array<Mesh> meshes;
    // Array<Curves> curves;
    std::vector<Mesh> meshes;
    std::vector<Curves> curves;
    Device *device;

    Scene() = default;
    ~Scene() = default;
};

YBI_NAMESPACE_END
