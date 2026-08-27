#include "scene/scene.h"
#include "core/arena.h"
#include "core/assert.h"
#include "core/vec2f.h"
#include "device/device.h"
#include <cstdint>
#include <cstring>
#include <memory>
#include <unordered_map>
#include <vector>

namespace ybi
{

Instance::Instance(AffineTransform transform, GeometryHandle handle)
    : transform(transform), handle(handle)
{
}

Mesh::Mesh(BufferHandle positions,
           BufferHandle normals,
           AttributeType normalType,
           BufferHandle uvs,
           BufferHandle uvIndices,
           AttributeType uvType,
           BufferHandle indices)
    : positions(positions), normals(normals), uvs(uvs), uvIndices(uvIndices), indices(indices),
      uvType(uvType), normalType(normalType)
{
}

BufferHandle Mesh::GetPositionsHandle() const
{
    return positions;
}

BufferHandle Mesh::GetNormalsHandle() const
{
    return normals;
}

BufferHandle Mesh::GetUVsHandle() const
{
    return uvs;
}

BufferHandle Mesh::GetUVIndicesHandle() const
{
    return uvIndices;
}

BufferHandle Mesh::GetIndicesHandle() const
{
    return indices;
}

BufferHandle Mesh::GetFaceMaterialIndicesHandle() const
{
    return faceMaterialIndices;
}

AttributeType Mesh::GetUVType() const
{
    return uvType;
}

AttributeType Mesh::GetNormalType() const
{
    return normalType;
}

uint32_t Mesh::GetMaterialIndex() const
{
    return materialIndex;
}

void Mesh::SetPositionsHandle(BufferHandle handle)
{
    positions = handle;
}

void Mesh::SetNormalsHandle(BufferHandle handle)
{
    normals = handle;
}

void Mesh::SetUVsHandle(BufferHandle handle)
{
    uvs = handle;
}

void Mesh::SetUVIndicesHandle(BufferHandle handle)
{
    uvIndices = handle;
}

void Mesh::SetIndicesHandle(BufferHandle handle)
{
    indices = handle;
}

void Mesh::SetFaceMaterialIndicesHandle(BufferHandle handle)
{
    faceMaterialIndices = handle;
}

void Mesh::SetMaterialIndex(uint32_t index)
{
    materialIndex = index;
}

bool Mesh::HasMaterial() const
{
    return materialIndex != ~0u || faceMaterialIndices.IsValid();
}

bool Mesh::HasSingleMaterial() const
{
    return materialIndex != ~0u;
}

// -----------------------------------------------------------------------------
// Attributes
// -----------------------------------------------------------------------------

size_t Attributes::AddAttribute(Array<Vec2f> &&array)
{
    // TODO: deduplication
    vector2s.emplace_back(std::move(array));

    size_t index = vector2s.size() - 1;
    return index;
}

// -----------------------------------------------------------------------------
// Data Container
// -----------------------------------------------------------------------------

Array<BVH> &DataContainer::GetBVHs(GeometryType type)
{
    return bvhHandles[(int)type];
}

const Array<Vec3f> &DataContainer::GetMeshPositions(const Mesh &mesh) const
{
    return meshPositions[mesh.GetPositionsHandle().index];
}

const Array<Vec3f> &DataContainer::GetMeshNormals(const Mesh &mesh) const
{
    return meshNormals[mesh.GetNormalsHandle().index];
}

const Array<int> &DataContainer::GetMeshIndices(const Mesh &mesh) const
{
    return meshTriangleIndices[mesh.GetIndicesHandle().index];
}

const Array<int> &DataContainer::GetMeshUVIndices(const Mesh &mesh) const
{
    return meshUVIndices[mesh.GetUVIndicesHandle().index];
}

// Scene
Scene::Scene(DataContainer *container) : container(container) {}

} // namespace ybi
