#pragma once

#include "util/array.h"
#include "util/float3x4.h"
#include <vector>

YBI_NAMESPACE_BEGIN

struct float3;
struct Device;

enum PrimitiveType : int
{
    PRIMITIVE_TYPE_TRIANGLES,
    PRIMITIVE_TYPE_CURVES,
    PRIMITIVE_TYPE_INSTANCES,
};

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

struct Primitive
{
};

struct CollectionRange
{
    Primitive *_begin;
    Primitive *_end;

    Primitive *begin() const
    {
        return _begin;
    }
    Primitive *end() const
    {
        return _end;
    }
};

struct ConstCollectionRange
{
    const Primitive *_begin;
    const Primitive *_end;

    const Primitive *begin() const
    {
        return _begin;
    }
    const Primitive *end() const
    {
        return _end;
    }
};

struct Scene
{
    BVH bvh;
    std::vector<Mesh> meshes;
    std::vector<Curves> curves;
    std::vector<Instances> instancesArray;

    Array<Primitive> primitives;
    Array<int> primitiveCollections;
    Device *device;

    Scene() = default;
    ~Scene() = default;

    int GetNumPrimitives(int collectionIndex) const;
    ConstCollectionRange GetPrimitivesInCollection(int collectionIndex) const;
};

YBI_NAMESPACE_END
