#pragma once

#include "core/affine_transform.h"
#include "core/array.h"
#include "core/hash.h"
#include "core/matrix4x4.h"
#include "device/device_memory.h"
#include "scene/bvh.h"
#include "scene/device_scene.h"
#include "scene/material.h"
#include "scene/openpbr_parameters.h"
#include <array>
#include <cstring>
#include <functional>
#include <limits>
#include <memory>
#include <string>
#include <unordered_map>
#include <variant>
#include <vector>

namespace ybi
{

struct Device;
class DeviceArena;

enum class PrimvarValueType
{
    Bool,
    Int,
    Float,
    ArrayVec2f,
    ArrayVec3f,
    ArrayVec4f,
};

using PrimvarValue = std::variant<bool, int, float, Array<Vec2f>, Array<Vec3f>, Array<Vec4f>>;

inline bool IsPrimvarArray(PrimvarValueType type)
{
    return type == PrimvarValueType::ArrayVec2f || type == PrimvarValueType::ArrayVec3f ||
           type == PrimvarValueType::ArrayVec4f;
}

struct Primvar
{
    AttributeId id;
    AttributeType type;
    PrimvarValueType valueType;
    size_t index;
    size_t indicesIndex = std::numeric_limits<size_t>::max();

    bool HasIndices() const;
};

struct MeshAttribute
{
    AttributeType type;
    PrimvarValue value;
    Array<int> indices;
};

struct Material
{
    std::string path;
    UdimTextureVector diffuseUdimPaths;
    uint32_t udimBase = 0;
    uint32_t udimTableOffset = 0;

    bool valid;
};

struct RectLight
{
    AffineTransform localToWorld;
    float width;
    float height;
    Vec3f radiance; // TODO: this can be a texture too
    uint32_t meshIndex = ~0u;
};

enum class GeometryType
{
    Mesh,
    Subdivision,
    Curve,
    Scene,
    Max,
};

inline const char *ConvertGeometryTypeToString(GeometryType type)
{
    switch (type)
    {
        case GeometryType::Mesh:
        {
            return "Mesh";
        }
        case GeometryType::Subdivision:
        {
            return "Subdivision";
        }
        case GeometryType::Curve:
        {
            return "Curve";
        }
        case GeometryType::Scene:
        {
            return "Scene";
        }
        default:
        {
            return "Invalid";
        }
    }
}

struct GeometryHandle
{
    GeometryType type;
    int index;
};

struct BufferHandle
{
    uint32_t index = ~0u;

    bool IsValid() const
    {
        return index != ~0u;
    }
};

struct Instance
{
    AffineTransform transform;
    GeometryHandle handle;

    Instance() = default;
    Instance(AffineTransform transform, GeometryHandle handle);
};

class Mesh
{
    BufferHandle positions;
    BufferHandle normals;
    BufferHandle uvs;
    BufferHandle uvIndices;
    BufferHandle indices;
    BufferHandle faceMaterialIndices;
    uint32_t materialIndex = ~0u;
    AttributeType uvType = AttributeType::None;
    AttributeType normalType = AttributeType::None;
    std::vector<Primvar> primvars;

public:
    Mesh() = default;
    Mesh(BufferHandle positions,
         BufferHandle normals,
         AttributeType normalType,
         BufferHandle uvs,
         BufferHandle uvIndices,
         AttributeType uvType,
         BufferHandle indices);

    BufferHandle GetPositionsHandle() const;
    BufferHandle GetNormalsHandle() const;
    BufferHandle GetUVsHandle() const;
    BufferHandle GetUVIndicesHandle() const;
    BufferHandle GetIndicesHandle() const;
    BufferHandle GetFaceMaterialIndicesHandle() const;
    AttributeType GetUVType() const;
    AttributeType GetNormalType() const;
    uint32_t GetMaterialIndex() const;
    void SetPositionsHandle(BufferHandle handle);
    void SetNormalsHandle(BufferHandle handle);
    void SetUVsHandle(BufferHandle handle);
    void SetUVIndicesHandle(BufferHandle handle);
    void SetIndicesHandle(BufferHandle handle);
    void SetFaceMaterialIndicesHandle(BufferHandle handle);
    void SetMaterialIndex(uint32_t index);
    void SetPrimvar(Primvar primvar);
    void SortPrimvars();
    bool PrimvarsAreSorted() const;
    const Primvar *GetPrimvar(AttributeId id) const;
    const std::vector<Primvar> &GetPrimvars() const;
    bool HasMaterial() const;
    bool HasSingleMaterial() const;
};

class Curve
{
};

class SubdivisionMesh
{
};

struct Camera
{
    Vec3f position;
    AffineTransform cameraFromWorld;
    Matrix4x4 clipFromCamera;
    Matrix4x4 rasterFromCamera;
    float verticalFovRadians;

    int width;
    int height;
};

struct Attributes
{
    std::vector<bool> bools;
    std::vector<int> ints;
    std::vector<float> floats;
    std::vector<Array<Vec2f>> vector2s;
    std::vector<Array<Vec3f>> vector3s;
    std::vector<Array<Vec4f>> vector4s;
    std::vector<Array<int>> indices;

    std::unordered_multimap<uint64_t, size_t> vector2HashToIndex;
    std::unordered_multimap<uint64_t, size_t> vector3HashToIndex;
    std::unordered_multimap<uint64_t, size_t> vector4HashToIndex;
    std::unordered_multimap<uint64_t, size_t> indicesHashToIndex;

    Primvar AddAttribute(AttributeId id,
                         AttributeType type,
                         PrimvarValue &&value,
                         Array<int> &&attributeIndices);

    template <typename T>
    const std::vector<Array<T>> &GetArrays() const
    {
        if constexpr (std::is_same_v<T, Vec2f>)
        {
            return vector2s;
        }
        else if constexpr (std::is_same_v<T, Vec3f>)
        {
            return vector3s;
        }
        else if constexpr (std::is_same_v<T, Vec4f>)
        {
            return vector4s;
        }
        else
        {
            static_assert("Unsupported type");
        }
    }

    const std::vector<Array<int>> &GetIndices() const;
    const Array<int> &GetIndices(size_t index) const;
};

using GeometryAttributeLoader = std::function<bool(const std::string &, MeshAttribute &)>;

struct PendingMesh
{
    size_t meshIndex;
    // std::vector<GeometryAttributeRequirement> requiredAttributes;
    GeometryAttributeLoader loadAttribute;
};

// NOTE: temporary context
struct SceneImportContext
{
    std::vector<PendingMesh> pendingMeshes;
};

// TODO: this is maybe getting messy
struct DataContainer
{
    Camera camera;
    std::vector<Mesh> meshes;
    std::vector<Array<Vec3f>> meshPositions;
    std::vector<Array<Vec3f>> meshNormals;
    std::vector<Array<Vec2f>> meshUVs;
    std::vector<Array<int>> meshUVIndices;
    std::vector<Array<int>> meshTriangleIndices;
    std::vector<Array<int>> meshFaceMaterialIndices;
    std::vector<Array<uint8_t>> compressedClusterBytes;
    std::vector<Material> materials;
    Curve *curves;
    SubdivisionMesh *subdivisionMeshes;

    std::vector<OpenPbrMaterial> openPbrMaterials;
    std::vector<MaterialHandle> materialHandles;

    std::vector<RectLight> rectLights;

    Attributes attributes;

    // BVHs
    std::array<Array<BVH>, int(GeometryType::Max)> bvhHandles;

    // Default handles
    BufferHandle defaultRectHandle;
    BVH defaultRectBVH;

    void Init(Device *device, DeviceArena *deviceArena);
    Array<BVH> &GetBVHs(GeometryType type);
    const Array<Vec3f> &GetMeshPositions(const Mesh &mesh) const;
    const Array<Vec3f> &GetMeshNormals(const Mesh &mesh) const;
    const Array<int> &GetMeshIndices(const Mesh &mesh) const;
    const Array<int> &GetMeshUVIndices(const Mesh &mesh) const;
};

struct Scene
{
    DataContainer *container;

    // NOTE: handle points to geometry in DataContainer. DataContainer may be shared
    // b/t multiple scenes
    Array<Instance> instances;
    // NOTE: these cannot be GeometryType::Scene
    // Array<GeometryHandle> geometryHandles;

    Scene() = default;
    Scene(DataContainer *container);
};

struct DataCache
{
    std::unordered_multimap<uint64_t, size_t> hashToIndex;

    DataCache() = default;
    ~DataCache() = default;

    template <typename T, typename GetFn>
    size_t GetOrAddArray(const Array<T> &array, size_t storedIndex, GetFn get)
    {
        const size_t numBytes = sizeof(T) * array.size();
        const uint8_t empty = 0;
        const uint8_t *bytes = numBytes ? reinterpret_cast<const uint8_t *>(array.data()) : &empty;
        uint64_t hash = MurmurHash64A(bytes, numBytes, 0);

        auto range = hashToIndex.equal_range(hash);
        for (auto it = range.first; it != range.second; ++it)
        {
            size_t index = it->second;
            const Array<T> &otherArray = get(index);
            if (otherArray.size() != array.size())
            {
                continue;
            }

            int result = numBytes ? memcmp(array.data(), otherArray.data(), numBytes) : 0;
            if (result == 0)
            {
                return index;
            }
        }

        hashToIndex.emplace(hash, storedIndex);
        return ~0ull;
    }
};

} // namespace ybi
