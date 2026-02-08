#include "scene.h"
#include "util/float3.h"

YBI_NAMESPACE_BEGIN

BVH::BVH() : flags(BVHFlags(0)) {}

Mesh::Mesh(Array<float3> &&pos, Array<int> &&idx)
    : positions(std::move(pos)), indices(std::move(idx))
{
}

Curves::Curves(Array<float3> &&positions, Array<float> &&widths, Array<int> &&curveVertexOffsets)
    : positions(std::move(positions)), widths(std::move(widths)),
      curveVertexOffsets(std::move(curveVertexOffsets))
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
    assert(count >= 4);
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

YBI_NAMESPACE_END
