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

template <typename T, size_t alignment>
static Array<T, alignment> CopyArray(const Array<T, alignment> &src)
{
    Array<T, alignment> out(src.size());
    if (src.size() > 0)
    {
        std::memcpy(out.data(), src.data(), src.size() * sizeof(T));
    }
    return out;
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

static bool BuildTriangulatedCornerRemap(const SubdivisionMesh &mesh,
                                         std::vector<int> *outCornerRemap,
                                         std::string *error)
{
    YBI_ASSERT(outCornerRemap);
    outCornerRemap->clear();

    size_t cornerOffset = 0;
    for (size_t faceIndex = 0; faceIndex < mesh.vertsPerFace.size(); ++faceIndex)
    {
        const int faceCount = mesh.vertsPerFace[faceIndex];
        if (faceCount < 3)
        {
            SetError(error,
                     "FlattenScenePoolToRootChildren: subdivision face has fewer than 3 vertices");
            return false;
        }
        if (cornerOffset + static_cast<size_t>(faceCount) > mesh.indices.size())
        {
            SetError(error,
                     "FlattenScenePoolToRootChildren: subdivision topology exceeds corner index range");
            return false;
        }

        for (int i = 1; i + 1 < faceCount; ++i)
        {
            outCornerRemap->push_back(static_cast<int>(cornerOffset + 0));
            outCornerRemap->push_back(static_cast<int>(cornerOffset + i));
            outCornerRemap->push_back(static_cast<int>(cornerOffset + i + 1));
        }
        cornerOffset += static_cast<size_t>(faceCount);
    }

    if (cornerOffset != mesh.indices.size())
    {
        SetError(error,
                 "FlattenScenePoolToRootChildren: subdivision topology corner count mismatch");
        return false;
    }
    return true;
}

static Attribute CopySubdivisionAttributeToMesh(const Attribute &srcAttr,
                                                size_t cornerCount,
                                                const std::vector<int> &cornerRemap)
{
    Array<uint8_t> data = CopyArray(srcAttr.data);
    if (srcAttr.interpolation == PrimvarInterpolation::FaceVarying &&
        srcAttr.indices.size() == cornerCount)
    {
        Array<int> indices(cornerRemap.size());
        for (size_t i = 0; i < cornerRemap.size(); ++i)
        {
            indices[i] = srcAttr.indices[static_cast<size_t>(cornerRemap[i])];
        }
        return Attribute(
            std::move(data), std::move(indices), srcAttr.type, srcAttr.interpolation, srcAttr.name);
    }

    if (srcAttr.indices.size() > 0)
    {
        Array<int> indices = CopyArray(srcAttr.indices);
        return Attribute(
            std::move(data), std::move(indices), srcAttr.type, srcAttr.interpolation, srcAttr.name);
    }

    return Attribute(std::move(data), srcAttr.type, srcAttr.interpolation, srcAttr.name);
}

static bool CopySubdivisionMeshToMesh(const SubdivisionMesh &src, Mesh *out, std::string *error)
{
    YBI_ASSERT(out);

    std::vector<int> cornerRemap;
    if (!BuildTriangulatedCornerRemap(src, &cornerRemap, error))
    {
        return false;
    }

    Array<int> triangulatedIndices(cornerRemap.size());
    for (size_t i = 0; i < cornerRemap.size(); ++i)
    {
        triangulatedIndices[i] = src.indices[static_cast<size_t>(cornerRemap[i])];
    }

    out->positions = CopyArray(src.vertices);
    out->indices = std::move(triangulatedIndices);
    out->attributes.clear();
    out->attributes.reserve(src.attributes.size());
    for (const Attribute &attr : src.attributes)
    {
        out->attributes.push_back(
            CopySubdivisionAttributeToMesh(attr, src.indices.size(), cornerRemap));
    }
    out->parentFromLocal = src.parentFromLocal;
    out->primPath = src.primPath;
    out->materialIndex = src.materialIndex;
    out->refIndex = UINT32_MAX;
    return true;
}

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
    dst->materials = src->materials;
    dst->lights = src->lights;

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

    auto getOrCreateLeafForSubdivision = [&](Scene *scene,
                                             size_t subdivIndex,
                                             uint32_t *outLeafIndex) -> bool {
        YBI_ASSERT(outLeafIndex);
        const PrimitiveKey key = {scene, 2u, static_cast<uint32_t>(subdivIndex)};
        auto found = leafSceneByPrimitive.find(key);
        if (found != leafSceneByPrimitive.end())
        {
            *outLeafIndex = found->second;
            return true;
        }

        std::unique_ptr<Scene> leaf = std::make_unique<Scene>();
        leaf->meshes.emplace_back();
        if (!CopySubdivisionMeshToMesh(scene->subdivisionMeshes[subdivIndex], &leaf->meshes.back(), error))
        {
            return false;
        }
        leaf->meshes.back().parentFromLocal = GetIdentityTransform();
        const uint32_t leafIndex = static_cast<uint32_t>(dst->scenes.size());
        dst->scenes.push_back(std::move(leaf));
        leafSceneByPrimitive.emplace(key, leafIndex);
        *outLeafIndex = leafIndex;
        return true;
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

        for (size_t subdivIndex = 0; subdivIndex < scene->subdivisionMeshes.size(); ++subdivIndex)
        {
            const float4x4 worldFromMesh = mul(
                entry.worldFromScene, ToFloat4x4(scene->subdivisionMeshes[subdivIndex].parentFromLocal));
            uint32_t leafIndex = 0u;
            if (!getOrCreateLeafForSubdivision(scene, subdivIndex, &leafIndex))
            {
                return false;
            }
            Scene *leafScene = dst->scenes[leafIndex].get();
            YBI_ASSERT(leafScene);
            const uint32_t rootChildIndex = static_cast<uint32_t>(dstRoot->childScenes.size());
            dstRoot->childScenes.push_back(leafScene);
            dstRoot->instances.emplace_back(ToFloat3x4(worldFromMesh), rootChildIndex);
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
