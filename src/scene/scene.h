#pragma once

#include "scene/attributes.h"
#include "scene/micropolygon_mesh.h"
#include "scene/subdivision_mesh.h"
#include "util/array.h"
#include "util/float2.h"
#include "util/float3.h"
#include "util/float3x4.h"
#include "util/float4x4.h"
#include "util/host_memory_arena.h"
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

YBI_NAMESPACE_BEGIN

struct float3;
struct Device;

enum PrimitiveType : int
{
    PRIMITIVE_TYPE_TRIANGLES,
    PRIMITIVE_TYPE_CURVES,
    PRIMITIVE_TYPE_INSTANCES,
    PRIMITIVE_TYPE_MAX,
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
    std::vector<Attribute> attributes;
    float3x4 parentFromLocal;
    uint32_t refIndex = UINT32_MAX;
    int materialIndex = -1;

    Mesh() = default;
    ~Mesh() = default;
    Mesh(Mesh &&other) = default;

    Mesh(Array<float3> &&pos, Array<int> &&idx);
    Mesh(Array<float3> &&pos, Array<int> &&idx, const float3x4 &parentFromLocal);
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
    float3x4 parentFromLocal;

    Curves() = default;
    ~Curves() = default;
    Curves(Curves &&other) = default;

    Curves(Array<float3> &&positions, Array<float> &&widths, Array<int> &&curveVertexOffsets);
    Curves(Array<float3> &&positions,
           Array<float> &&widths,
           Array<int> &&curveVertexOffsets,
           const float3x4 &parentFromLocal);
    size_t GetNumVertices() const;
    size_t GetNumCurves() const;
    size_t GetNumSegments() const;
    int GetCurveKeyStart(size_t curveIndex) const;
    int GetCurveNumSegments(size_t curveIndex) const;
    void GetCurveRange(uint32_t index, uint32_t &start, uint32_t &count) const;

    const Array<float3> &GetVertices() const;
    const Array<float> &GetWidths() const;
};

struct Instance
{
    float3x4 parentFromLocal;
    uint32_t childSceneIndex;

    Instance() = default;
    ~Instance() = default;
    Instance(Instance &&other) = default;

    Instance(const float3x4 &parentFromLocal, uint32_t childSceneIndex);
};

struct Primitive
{
    PrimitiveType primitiveType;
    int index;
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

struct Camera
{
    int viewportWidth = 1920;
    int viewportHeight = 1080;
    float4x4 cameraFromWorld = float4x4::Identity();
    float4x4 clipFromCamera = float4x4::Identity();
    float3 worldPosition = make_float3(0.0f, 0.0f, 0.0f);
    float3 forward = make_float3(0.0f, 0.0f, -1.0f);
    float verticalFovDegrees = 45.0f;
    float nearPlane = 1.0f;
    bool hasValidCamera = false;
    std::string path;
};

struct Scene
{
    using BVHHandle = uint64_t;

    BVH bvh;
    BVHHandle bvhHandle = 0;
    std::vector<Mesh> meshes;
    std::vector<Curves> curves;
    std::vector<Instance> instances;
    std::vector<Scene *> childScenes;
    std::vector<MicropolygonMesh> micropolygonMeshes;

    // TODO: these are tessellated and displaced to either Mesh or Micropolygon Mesh
    std::vector<SubdivisionMesh> subdivisionMeshes;

    Array<Primitive> primitives;
    Array<int> primitiveCollections;

    // HostMemoryArena arena;
    Array<Attribute> attributes;
    Device *device;

    Scene() = default;
    ~Scene() = default;
    Scene(const Scene &) = delete;
    Scene &operator=(const Scene &) = delete;
    Scene(Scene &&other) noexcept;
    Scene &operator=(Scene &&other) noexcept = delete;
};

struct ScenePool
{
    struct MaterialTextureInput
    {
        std::string inputName;
        std::string texturePath;
        std::string swizzle;
    };

    struct MaterialInfo
    {
        std::string materialPath;
        std::vector<MaterialTextureInput> textureInputs;
        std::string ntcDiffuseFile;
        std::string ntcDiffuseTextureName;
    };

    std::vector<std::unique_ptr<Scene>> scenes;
    uint32_t rootSceneIndex = 0;
    Camera camera;
    std::vector<MaterialInfo> materials;
};

struct SceneMeshUploadRef
{
    const Mesh *mesh = nullptr;
    uint32_t refIndex = 0;
};

bool FlattenScenePoolToRootChildren(ScenePool *src, ScenePool *dst, std::string *error = nullptr);
void CollectScenePoolMeshUploadRefs(ScenePool *scenePool, std::vector<SceneMeshUploadRef> *out);

YBI_NAMESPACE_END
