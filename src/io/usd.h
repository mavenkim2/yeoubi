#pragma once

#include "core/affine_transform.h"
#include "core/array.h"
#include "core/concurrent_hash_map.h"
#include "core/matrix4x4.h"
#include "scene/scene.h"
#include <atomic>
#include <cstdint>
#include <memory>
#include <pxr/base/tf/token.h>
#include <pxr/base/vt/array.h>
#include <pxr/pxr.h>
#include <pxr/usd/usd/prim.h>
#include <string>
#include <vector>

PXR_NAMESPACE_OPEN_SCOPE
class UsdGeomMesh;
class UsdGeomSubset;
class UsdGeomXformCache;
class UsdShadeMaterial;
PXR_NAMESPACE_CLOSE_SCOPE

namespace ybi
{

enum class AttributeType;
class Arena;
struct DataContainer;
struct Scene;
struct USDPrimTraversalData;

class SceneLoader
{
public:
    static std::unique_ptr<SceneLoader> CreateSceneLoader(const std::string &filePath);
    virtual Scene *LoadScene(Arena *arena, DataContainer *container) = 0;
};

template <typename T>
struct InternedMeshBuffer
{
    uint64_t hash = 0;
    Array<T> buffer;
};

template <typename T>
struct MeshBufferInternState
{
    ConcurrentHashMap map;
    std::atomic<uint32_t> bufferOffset;
    std::vector<InternedMeshBuffer<T>> buffers;
    Array<size_t> indices;
    Array<uint32_t> oldToNew;

    MeshBufferInternState(Arena *arena, size_t numMeshes);
    size_t Intern(Arena *arena, Array<T> &&array, size_t bufferIndex);
    void MoveSortedBuffersTo(std::vector<Array<T>> *dst);
    uint32_t ResolveIndex(size_t index) const;
};

struct MeshBufferInternStates
{
    MeshBufferInternState<Vec3f> positions;
    MeshBufferInternState<Vec3f> normals;
    MeshBufferInternState<Vec2f> uvs;
    MeshBufferInternState<int> indices;
    MeshBufferInternState<int> faceMaterialIndices;

    MeshBufferInternStates(Arena *arena, size_t numMeshes);
};

class USDSceneLoader : public SceneLoader
{
    std::string filePath;
    std::vector<pxr::TfToken> allowedPurposes;

    AffineTransform ConvertAffineTransform(const pxr::GfMatrix4d &m) const;
    Matrix4x4 ConvertMatrix4x4(const pxr::GfMatrix4d &m) const;
    void ConvertMeshIndicesToTriangles(const pxr::VtIntArray &faceIndices,
                                       const pxr::VtIntArray &faceCounts,
                                       Array<int> &triangles) const;
    void ConvertFaceVaryingUVsToTriangles(const pxr::VtVec2fArray &values,
                                          const pxr::VtIntArray &faceCounts,
                                          Array<Vec2f> &triangles) const;
    void ConvertFaceVaryingNormalsToTriangles(const pxr::VtVec3fArray &values,
                                              const pxr::VtIntArray &faceCounts,
                                              Array<Vec3f> &triangles) const;
    AttributeType ProcessPrimvars(const pxr::UsdGeomMesh &mesh,
                                  pxr::VtArray<pxr::GfVec2f> *outUVs,
                                  pxr::UsdTimeCode timeCode) const;
    AttributeType ProcessNormals(const pxr::UsdGeomMesh &mesh,
                                 pxr::VtVec3fArray *outNormals,
                                 pxr::UsdTimeCode timeCode) const;
    void ProcessPolygonMesh(Scene *scene,
                            pxr::UsdPrim sceneRoot,
                            size_t meshOffset,
                            size_t meshIndex,
                            const pxr::UsdPrim prim,
                            pxr::UsdGeomXformCache &cache,
                            pxr::UsdTimeCode timeCode = pxr::UsdTimeCode::Default()) const;
    void ProcessPolygonMeshConcurrent(Arena *arena,
                                      Scene *scene,
                                      pxr::UsdPrim sceneRoot,
                                      size_t meshOffset,
                                      size_t meshIndex,
                                      const pxr::UsdPrim prim,
                                      MeshBufferInternStates *internStates,
                                      pxr::UsdGeomXformCache &cache,
                                      pxr::UsdTimeCode timeCode);
    void GetFaceMaterialIndices(Array<int> &faceMaterialIndices,
                                const pxr::VtIntArray &faceCounts,
                                const std::vector<pxr::UsdGeomSubset> &subsets,
                                std::unordered_map<std::string, size_t> &materialMap,
                                std::vector<pxr::UsdShadeMaterial> &materials,
                                pxr::UsdTimeCode timeCode) const;
    void TraversePrim(USDPrimTraversalData *data,
                      pxr::UsdPrim root,
                      bool skipPrototypePrims = true) const;
    bool IsPurposeAllowed(pxr::UsdPrim prim) const;

public:
    USDSceneLoader(std::string filePath);
    ~USDSceneLoader() = default;
    Scene *LoadScene(Arena *arena, DataContainer *container) override;
};

} // namespace ybi
