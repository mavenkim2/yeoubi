#pragma once

#include "util/array.h"
#include "util/float3x4.h"
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
private:
    Array<float3> positions;
    Array<float> widths;
    // num_offsets = num_curves + 1
    Array<int> curveVertexOffsets;
    int curveFlags;

public:
    Curves() = default;
    ~Curves() = default;
    Curves(Curves &&other) = default;

    Curves(Array<float3> &&positions, Array<float> &&widths, Array<int> &&curveVertexOffsets);
    size_t GetNumVertices() const;
    size_t GetNumCurves() const;
    size_t GetNumSegments() const;
    int GetCurveKeyStart(size_t curveIndex) const;
    int GetCurveNumSegments(size_t curveIndex) const;
    void GetCurveRange(uint32_t index, uint32_t &start, uint32_t &count) const;

    const Array<float3> &GetVertices() const;
    const Array<float> &GetWidths() const;
};

struct Instances
{
private:
    Array<float3x4> affineTransforms;
    Array<int> objectIDs;

public:
    Instances() = default;
    ~Instances() = default;
    Instances(Instances &&other) = default;

    Instances(Array<float3x4> &&affineTransforms, Array<int> &&objectIDs);
};

struct Scene
{
    BVH bvh;
    std::vector<Mesh> meshes;
    std::vector<Curves> curves;
    std::vector<Instances> instancesArray;
    Device *device;

    Scene() = default;
    ~Scene() = default;
};

YBI_NAMESPACE_END
