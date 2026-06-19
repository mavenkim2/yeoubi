#pragma once

#include "core/affine_transform.h"
#include "core/array.h"
#include "scene/bvh.h"
#include <array>
#include <memory>
#include <vector>

namespace ybi
{

struct Device;
class DeviceArena;

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
    Array<int> indices;

public:
    Mesh() = default;
    Mesh(Array<Vec3f> &&positions);
    Mesh(Array<Vec3f> &&positions, Array<int> &&indices);
    Mesh(Mesh &&other) = default;
    Mesh &operator=(Mesh &&other) = default;

    size_t GetNumVertices() const;
    const Array<Vec3f> &GetPositions() const;
    const Array<int> &GetIndices() const;
};

class Curve
{
};

class SubdivisionMesh
{
};

struct DataContainer
{
    std::vector<Mesh> meshes;
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
    void
    Update(Device *device, DeviceArena *deviceArena, DataContainer *container, Scene *rootScene);
};

struct RenderSystem
{
    Scene *rootScene;
    DataContainer *dataContainer;

    Device *device;
    std::unique_ptr<DeviceArena> deviceArena;

    // Data managers
    GeometryManager geometryManager;

    RenderSystem(Device *device, DataContainer *container, Scene *rootScene);
    void Update();
};

} // namespace ybi
