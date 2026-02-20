#include "scene.h"
#include "util/assert.h"
#include "util/float3.h"

YBI_NAMESPACE_BEGIN

BVH::BVH() : flags(BVHFlags(0)) {}

static float3x4 GetIdentityTransform()
{
    return float3x4(1.0f,
                    0.0f,
                    0.0f,
                    0.0f,
                    0.0f,
                    1.0f,
                    0.0f,
                    0.0f,
                    0.0f,
                    0.0f,
                    1.0f,
                    0.0f);
}

Mesh::Mesh(Array<float3> &&pos, Array<int> &&idx)
    : positions(std::move(pos)), indices(std::move(idx)), parentFromLocal(GetIdentityTransform())
{
}

Mesh::Mesh(Array<float3> &&pos, Array<int> &&idx, const float3x4 &parentFromLocal)
    : positions(std::move(pos)), indices(std::move(idx)), parentFromLocal(parentFromLocal)
{
}

Curves::Curves(Array<float3> &&positions, Array<float> &&widths, Array<int> &&curveVertexOffsets)
    : positions(std::move(positions)), widths(std::move(widths)),
      curveVertexOffsets(std::move(curveVertexOffsets)), parentFromLocal(GetIdentityTransform())
{
}

Curves::Curves(Array<float3> &&positions,
               Array<float> &&widths,
               Array<int> &&curveVertexOffsets,
               const float3x4 &parentFromLocal)
    : positions(std::move(positions)), widths(std::move(widths)),
      curveVertexOffsets(std::move(curveVertexOffsets)), parentFromLocal(parentFromLocal)
{
}

size_t Curves::GetNumVertices() const
{
    return positions.size();
}

size_t Curves::GetNumCurves() const
{
    return curveVertexOffsets.size();
}

size_t Curves::GetNumSegments() const
{
    return positions.size() - 3 * GetNumCurves();
}

int Curves::GetCurveKeyStart(size_t curveIndex) const
{
    return curveVertexOffsets[curveIndex];
}

int Curves::GetCurveNumSegments(size_t curveIndex) const
{
    uint32_t start = curveVertexOffsets[curveIndex];
    uint32_t count = curveIndex == curveVertexOffsets.size() - 1
                         ? positions.size()
                         : curveVertexOffsets[curveIndex + 1];
    count -= start;
    YBI_ASSERT(count >= 4);
    return count - 3;
}

void Curves::GetCurveRange(uint32_t index, uint32_t &start, uint32_t &count) const
{
    start = curveVertexOffsets[index];
    count =
        index == curveVertexOffsets.size() - 1 ? positions.size() : curveVertexOffsets[index + 1];
    count -= start;
}

const Array<float3> &Curves::GetVertices() const
{
    return positions;
}

const Array<float> &Curves::GetWidths() const
{
    return widths;
}

Instance::Instance(const float3x4 &parentFromLocal, uint32_t childSceneIndex)
    : parentFromLocal(parentFromLocal), childSceneIndex(childSceneIndex)
{
}

Scene::Scene(Scene &&other) noexcept
    : bvh(other.bvh), bvhHandle(other.bvhHandle),
      meshes(std::move(other.meshes)),
      curves(std::move(other.curves)),
      instances(std::move(other.instances)),
      childScenes(std::move(other.childScenes)),
      micropolygonMeshes(std::move(other.micropolygonMeshes)),
      subdivisionMeshes(std::move(other.subdivisionMeshes)),
      primitives(std::move(other.primitives)),
      primitiveCollections(std::move(other.primitiveCollections)),
      attributes(std::move(other.attributes)), device(other.device)
{
}

int Scene::GetNumPrimitives(int collectionIndex) const
{
    return primitiveCollections[collectionIndex + 1] - primitiveCollections[collectionIndex];
}

ConstCollectionRange Scene::GetPrimitivesInCollection(int collectionIndex) const
{
    if (collectionIndex < 0 || collectionIndex >= (int)primitiveCollections.size() - 1)
    {
        return {nullptr, nullptr};
    }

    int startOffset = primitiveCollections[collectionIndex];
    int endOffset = primitiveCollections[collectionIndex + 1];

    return {&primitives[startOffset], &primitives[endOffset]};
}

YBI_NAMESPACE_END
