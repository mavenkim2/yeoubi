#include "scene.h"
#include "util/assert.h"
#include "util/float3.h"
#include "util/float4x4.h"

#include <string>
#include <unordered_map>
#include <vector>

YBI_NAMESPACE_BEGIN

BVH::BVH() : flags(BVHFlags(0)) {}

static float3x4 GetIdentityTransform()
{
    return float3x4(1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f);
}

static void SetError(std::string *error, const std::string &message)
{
    if (error)
    {
        *error = message;
    }
}

static float4x4 ToFloat4x4(const float3x4 &m)
{
    return float4x4(m.m[0][0],
                    m.m[0][1],
                    m.m[0][2],
                    m.m[0][3],
                    m.m[1][0],
                    m.m[1][1],
                    m.m[1][2],
                    m.m[1][3],
                    m.m[2][0],
                    m.m[2][1],
                    m.m[2][2],
                    m.m[2][3],
                    0.0f,
                    0.0f,
                    0.0f,
                    1.0f);
}

static float3x4 ToFloat3x4(const float4x4 &m)
{
    return float3x4(m.m[0][0],
                    m.m[0][1],
                    m.m[0][2],
                    m.m[0][3],
                    m.m[1][0],
                    m.m[1][1],
                    m.m[1][2],
                    m.m[1][3],
                    m.m[2][0],
                    m.m[2][1],
                    m.m[2][2],
                    m.m[2][3]);
}

struct PrimitiveKey
{
    const Scene *scene = nullptr;
    uint32_t type = 0;
    uint32_t index = 0;

    bool operator==(const PrimitiveKey &other) const
    {
        return scene == other.scene && type == other.type && index == other.index;
    }
};

struct PrimitiveKeyHash
{
    size_t operator()(const PrimitiveKey &k) const
    {
        const size_t a = std::hash<const void *>{}(k.scene);
        const size_t b = std::hash<uint32_t>{}(k.type);
        const size_t c = std::hash<uint32_t>{}(k.index);
        return a ^ (b + 0x9e3779b97f4a7c15ull + (a << 6) + (a >> 2)) ^
               (c + 0x9e3779b97f4a7c15ull + (b << 6) + (b >> 2));
    }
};

struct FlattenTraversalEntry
{
    Scene *scene = nullptr;
    float4x4 worldFromScene = float4x4::Identity();
};

bool FlattenScenePoolToRootChildren(ScenePool *src, ScenePool *dst, std::string *error)
{
    if (!src || !dst)
    {
        SetError(error, "FlattenScenePoolToRootChildren: null ScenePool");
        return false;
    }
    if (src == dst)
    {
        SetError(error, "FlattenScenePoolToRootChildren: src and dst must differ");
        return false;
    }
    if (src->scenes.empty())
    {
        SetError(error, "FlattenScenePoolToRootChildren: source has no scenes");
        return false;
    }
    if (src->rootSceneIndex >= src->scenes.size())
    {
        SetError(error, "FlattenScenePoolToRootChildren: source rootSceneIndex out of range");
        return false;
    }

    Scene *srcRoot = src->scenes[src->rootSceneIndex].get();
    if (!srcRoot)
    {
        SetError(error, "FlattenScenePoolToRootChildren: source root scene is null");
        return false;
    }

    dst->scenes.clear();
    dst->rootSceneIndex = 0;
    dst->camera = src->camera;

    std::unique_ptr<Scene> dstRoot = std::make_unique<Scene>();
    dstRoot->meshes = std::move(srcRoot->meshes);
    dstRoot->curves = std::move(srcRoot->curves);
    dstRoot->subdivisionMeshes = std::move(srcRoot->subdivisionMeshes);

    std::unordered_map<PrimitiveKey, uint32_t, PrimitiveKeyHash> leafSceneByPrimitive;
    leafSceneByPrimitive.reserve(src->scenes.size() * 4);

    auto getOrCreateLeafForMesh = [&](Scene *scene, size_t meshIndex) -> uint32_t {
        const PrimitiveKey key = {scene, 0u, static_cast<uint32_t>(meshIndex)};
        auto found = leafSceneByPrimitive.find(key);
        if (found != leafSceneByPrimitive.end())
        {
            return found->second;
        }

        std::unique_ptr<Scene> leaf = std::make_unique<Scene>();
        leaf->meshes.emplace_back(std::move(scene->meshes[meshIndex]));
        leaf->meshes.back().parentFromLocal = GetIdentityTransform();
        const uint32_t leafIndex = static_cast<uint32_t>(dst->scenes.size());
        dst->scenes.push_back(std::move(leaf));
        leafSceneByPrimitive.emplace(key, leafIndex);
        return leafIndex;
    };

    auto getOrCreateLeafForCurves = [&](Scene *scene, size_t curveIndex) -> uint32_t {
        const PrimitiveKey key = {scene, 1u, static_cast<uint32_t>(curveIndex)};
        auto found = leafSceneByPrimitive.find(key);
        if (found != leafSceneByPrimitive.end())
        {
            return found->second;
        }

        std::unique_ptr<Scene> leaf = std::make_unique<Scene>();
        leaf->curves.emplace_back(std::move(scene->curves[curveIndex]));
        leaf->curves.back().parentFromLocal = GetIdentityTransform();
        const uint32_t leafIndex = static_cast<uint32_t>(dst->scenes.size());
        dst->scenes.push_back(std::move(leaf));
        leafSceneByPrimitive.emplace(key, leafIndex);
        return leafIndex;
    };

    std::vector<FlattenTraversalEntry> stack;
    stack.push_back({srcRoot, float4x4::Identity()});

    while (!stack.empty())
    {
        const FlattenTraversalEntry entry = stack.back();
        stack.pop_back();

        Scene *scene = entry.scene;
        YBI_ASSERT(scene);

        for (const Instance &instance : scene->instances)
        {
            if (instance.childSceneIndex >= scene->childScenes.size())
            {
                SetError(error,
                         "FlattenScenePoolToRootChildren: instance childSceneIndex out of range");
                return false;
            }
            Scene *child = scene->childScenes[instance.childSceneIndex];
            if (!child)
            {
                SetError(error, "FlattenScenePoolToRootChildren: null child scene");
                return false;
            }
            const float4x4 worldFromChild =
                mul(entry.worldFromScene, ToFloat4x4(instance.parentFromLocal));
            stack.push_back({child, worldFromChild});
        }

        if (scene == srcRoot)
        {
            continue;
        }

        for (size_t meshIndex = 0; meshIndex < scene->meshes.size(); meshIndex++)
        {
            const float4x4 worldFromMesh =
                mul(entry.worldFromScene, ToFloat4x4(scene->meshes[meshIndex].parentFromLocal));
            const uint32_t leafIndex = getOrCreateLeafForMesh(scene, meshIndex);
            Scene *leafScene = dst->scenes[leafIndex].get();
            YBI_ASSERT(leafScene);
            const uint32_t rootChildIndex = static_cast<uint32_t>(dstRoot->childScenes.size());
            dstRoot->childScenes.push_back(leafScene);
            dstRoot->instances.emplace_back(ToFloat3x4(worldFromMesh), rootChildIndex);
        }

        for (size_t curveIndex = 0; curveIndex < scene->curves.size(); curveIndex++)
        {
            const float4x4 worldFromCurve =
                mul(entry.worldFromScene, ToFloat4x4(scene->curves[curveIndex].parentFromLocal));
            const uint32_t leafIndex = getOrCreateLeafForCurves(scene, curveIndex);
            Scene *leafScene = dst->scenes[leafIndex].get();
            YBI_ASSERT(leafScene);
            const uint32_t rootChildIndex = static_cast<uint32_t>(dstRoot->childScenes.size());
            dstRoot->childScenes.push_back(leafScene);
            dstRoot->instances.emplace_back(ToFloat3x4(worldFromCurve), rootChildIndex);
        }
    }

    dst->scenes.push_back(std::move(dstRoot));
    dst->rootSceneIndex = static_cast<uint32_t>(dst->scenes.size() - 1u);

    src->scenes.clear();
    src->rootSceneIndex = 0;
    return true;
}

void CollectScenePoolMeshUploadRefs(ScenePool *scenePool, std::vector<SceneMeshUploadRef> *out)
{
    YBI_ASSERT(scenePool);
    YBI_ASSERT(out);
    out->clear();

    uint32_t nextRefIndex = 0;
    for (std::unique_ptr<Scene> &scenePtr : scenePool->scenes)
    {
        Scene *scene = scenePtr.get();
        if (!scene)
        {
            continue;
        }
        for (Mesh &mesh : scene->meshes)
        {
            mesh.refIndex = nextRefIndex;
            out->push_back({&mesh, nextRefIndex});
            nextRefIndex++;
        }
    }
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
    : bvh(other.bvh), bvhHandle(other.bvhHandle), meshes(std::move(other.meshes)),
      curves(std::move(other.curves)), instances(std::move(other.instances)),
      childScenes(std::move(other.childScenes)),
      micropolygonMeshes(std::move(other.micropolygonMeshes)),
      subdivisionMeshes(std::move(other.subdivisionMeshes)),
      primitives(std::move(other.primitives)),
      primitiveCollections(std::move(other.primitiveCollections)),
      attributes(std::move(other.attributes)), device(other.device)
{
}

YBI_NAMESPACE_END
