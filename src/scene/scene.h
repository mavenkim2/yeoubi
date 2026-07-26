#pragma once

#include "core/affine_transform.h"
#include "core/array.h"
#include "core/hash.h"
#include "core/matrix4x4.h"
#include "device/device_memory.h"
#include "scene/bvh.h"
#include "scene/device_scene.h"
#include "scene/texture.h"
#include <array>
#include <cstring>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

namespace ybi
{

struct Device;
class DeviceArena;

struct Material
{
    std::string path;
    std::vector<std::pair<std::string, uint32_t>> diffuseUdimPaths;
    uint32_t udimBase = 0;
    uint32_t udimTableOffset = 0;
};

enum class GeometryType
{
    Mesh,
    Subdivison,
    Curve,
    Scene,
    Max,
};

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

    // BVHs
    std::array<Array<BVH>, int(GeometryType::Max)> bvhHandles;

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
