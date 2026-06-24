#pragma once

#include "core/affine_transform.h"
#include "core/array.h"
#include "core/matrix4x4.h"
#include "pxr/pxr.h"
#include <memory>
#include <pxr/base/vt/array.h>
#include <pxr/usd/usd/prim.h>
#include <string>

PXR_NAMESPACE_OPEN_SCOPE
class UsdGeomXformCache;
PXR_NAMESPACE_CLOSE_SCOPE

namespace ybi
{

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

class USDSceneLoader : public SceneLoader
{
    std::string filePath;

    AffineTransform ConvertAffineTransform(const pxr::GfMatrix4d &m) const;
    Matrix4x4 ConvertMatrix4x4(const pxr::GfMatrix4d &m) const;
    void ConvertMeshToTriangles(const pxr::VtIntArray &faceIndices,
                                const pxr::VtIntArray &faceCounts,
                                Array<int> &triangles) const;
    void ProcessPolygonMesh(Scene *scene,
                            size_t meshOffset,
                            size_t meshIndex,
                            const pxr::UsdPrim prim,
                            pxr::UsdGeomXformCache &cache,
                            pxr::UsdTimeCode timeCode = pxr::UsdTimeCode::Default()) const;
    void TraversePrim(USDPrimTraversalData *data,
                      pxr::UsdPrim root,
                      bool skipPrototypePrims = true) const;

public:
    USDSceneLoader(std::string filePath);
    ~USDSceneLoader() = default;
    Scene *LoadScene(Arena *arena, DataContainer *container) override;
};

} // namespace ybi
