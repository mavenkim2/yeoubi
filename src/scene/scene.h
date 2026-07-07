#pragma once

#include "core/affine_transform.h"
#include "core/array.h"
#include "core/matrix4x4.h"
#include "device/device_memory.h"
#include "scene/bvh.h"
#include "scene/device_scene.h"
#include "scene/texture.h"
#include <array>
#include <memory>
#include <string>
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

struct Instance
{
    AffineTransform transform;
    GeometryHandle handle;

    Instance() = default;
    Instance(AffineTransform transform, GeometryHandle handle);
};

class Mesh
{
    Array<Vec3f> positions;
    Array<Vec2f> uvs;
    Array<int> indices;
    uint32_t materialIndex = ~0u;
    AttributeType uvType = AttributeType::None;

public:
    Mesh() = default;
    Mesh(Array<Vec3f> &&positions, Array<Vec2f> &&uvs, AttributeType uvType, Array<int> &&indices);
    Mesh(Mesh &&other) = default;
    Mesh &operator=(Mesh &&other) = default;

    size_t GetNumVertices() const;
    size_t GetNumUVs() const;
    AttributeType GetUVType() const;
    uint32_t GetMaterialIndex() const;
    const Array<Vec3f> &GetPositions() const;
    const Array<Vec2f> &GetUVs() const;
    const Array<int> &GetIndices() const;
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

    int width;
    int height;
};

// TODO: this is maybe getting messy
struct DataContainer
{
    Camera camera;
    std::vector<Mesh> meshes;
    std::vector<Material> materials;
    Curve *curves;
    SubdivisionMesh *subdivisionMeshes;

    // BVHs
    std::array<Array<BVH>, int(GeometryType::Max)> bvhHandles;

    Array<BVH> &GetBVHs(GeometryType type);
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

struct GeometryManager
{
    void Update(Device *device,
                DeviceArena *deviceArena,
                DataContainer *container,
                Scene *rootScene,
                DeviceScene *deviceScene);
};

struct RenderSystem
{
    Scene *rootScene;
    DataContainer *dataContainer;

    DeviceScene deviceScene;

    Device *device;
    std::unique_ptr<DeviceArena> deviceArena;

    // Data managers
    GeometryManager geometryManager;

    RenderSystem(Device *device, DataContainer *container, Scene *rootScene);
    void Update();
};

} // namespace ybi
