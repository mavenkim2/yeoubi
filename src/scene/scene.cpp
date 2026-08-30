#include "scene/scene.h"
#include "core/arena.h"
#include "core/assert.h"
#include "core/vec2f.h"
#include "device/device.h"
#include <cstdint>
#include <cstring>
#include <memory>
#include <type_traits>
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

void Mesh::SetPrimvar(Primvar primvar)
{
    for (Primvar &existing : primvars)
    {
        if (existing.id == primvar.id)
        {
            existing = std::move(primvar);
            return;
        }
    }

    primvars.emplace_back(std::move(primvar));
}

const Primvar *Mesh::GetPrimvar(AttributeId id) const
{
    for (const Primvar &primvar : primvars)
    {
        if (primvar.id == id)
        {
            return &primvar;
        }
    }

    return nullptr;
}

const std::vector<Primvar> &Mesh::GetPrimvars() const
{
    return primvars;
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

bool Primvar::HasIndices() const
{
    return indicesIndex != std::numeric_limits<size_t>::max();
}

Primvar Attributes::AddAttribute(AttributeId id,
                                 AttributeType type,
                                 PrimvarValue &&value,
                                 Array<int> &&attributeIndices)
{
    // TODO: deduplication
    Primvar primvar{id, type, PrimvarValueType::Bool, 0};

    std::visit(
        [&](auto &&attribute) {
            using T = std::decay_t<decltype(attribute)>;
            if constexpr (std::is_same_v<T, bool>)
            {
                primvar.valueType = PrimvarValueType::Bool;
                bools.push_back(attribute);
                primvar.index = bools.size() - 1;
            }
            else if constexpr (std::is_same_v<T, int>)
            {
                primvar.valueType = PrimvarValueType::Int;
                ints.push_back(attribute);
                primvar.index = ints.size() - 1;
            }
            else if constexpr (std::is_same_v<T, float>)
            {
                primvar.valueType = PrimvarValueType::Float;
                floats.push_back(attribute);
                primvar.index = floats.size() - 1;
            }
            else if constexpr (std::is_same_v<T, Array<Vec2f>>)
            {
                primvar.valueType = PrimvarValueType::ArrayVec2f;
                vector2s.emplace_back(std::move(attribute));
                primvar.index = vector2s.size() - 1;
            }
            else if constexpr (std::is_same_v<T, Array<Vec3f>>)
            {
                primvar.valueType = PrimvarValueType::ArrayVec3f;
                vector3s.emplace_back(std::move(attribute));
                primvar.index = vector3s.size() - 1;
            }
            else if constexpr (std::is_same_v<T, Array<Vec4f>>)
            {
                primvar.valueType = PrimvarValueType::ArrayVec4f;
                vector4s.emplace_back(std::move(attribute));
                primvar.index = vector4s.size() - 1;
            }
        },
        std::move(value));

    if (attributeIndices.size())
    {
        indices.emplace_back(std::move(attributeIndices));
        primvar.indicesIndex = indices.size() - 1;
    }

    return primvar;
}

const std::vector<Array<int>> &Attributes::GetIndices() const
{
    return indices;
}

const Array<int> &Attributes::GetIndices(size_t index) const
{
    return indices[index];
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
