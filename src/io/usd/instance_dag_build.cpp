#include "io/usd/instance_dag_build.h"

#include "util/assert.h"
#include "util/float3.h"
#include "util/float3x4.h"
#include "util/float4.h"

#include <algorithm>
#include <cstdint>
#include <memory>
#include <pxr/base/gf/matrix4d.h>
#include <pxr/usd/sdf/path.h>
#include <pxr/usd/usd/prim.h>
#include <pxr/usd/usdGeom/basisCurves.h>
#include <pxr/usd/usdGeom/mesh.h>
#include <pxr/usd/usdGeom/pointInstancer.h>
#include <pxr/usd/usdGeom/xformCache.h>
#include <pxr/usd/usdGeom/xformable.h>
#include <string>
#include <unordered_map>
#include <vector>

YBI_NAMESPACE_BEGIN

namespace
{

struct USDPrimLists
{
    std::string ownerPath;
    std::vector<pxr::UsdPrim> meshes;
    std::vector<pxr::UsdPrim> curves;
    std::vector<pxr::UsdPrim> instances;
    std::vector<pxr::UsdPrim> pointInstancers;
};

struct SceneInstance
{
    float3x4 worldFromLocal;
    uint32_t childSceneIndex;
};

struct SceneMesh
{
    std::string path;
    float3x4 worldFromLocal;
};

struct SceneCurve
{
    std::string path;
    float3x4 worldFromLocal;
};

struct BuildScene
{
    std::string path;
    std::vector<BuildScene *> childScenes;
    std::vector<SceneMesh> meshes;
    std::vector<SceneCurve> curves;
    std::vector<SceneInstance> instances;
};

struct BuildResult
{
    std::unique_ptr<BuildScene> rootScene;
    std::vector<std::unique_ptr<BuildScene>> prototypeScenes;
};

static void SetError(std::string *error, const std::string &message)
{
    if (error)
    {
        *error = message;
    }
}

static float3x4 ConvertAffineTransform(const pxr::GfMatrix4d &m)
{
    const pxr::GfMatrix4d t = m.GetTranspose();
    const pxr::GfVec4d r0 = t.GetRow(0);
    const pxr::GfVec4d r1 = t.GetRow(1);
    const pxr::GfVec4d r2 = t.GetRow(2);

    return float3x4(make_float4((float)r0[0], (float)r0[1], (float)r0[2], (float)r0[3]),
                    make_float4((float)r1[0], (float)r1[1], (float)r1[2], (float)r1[3]),
                    make_float4((float)r2[0], (float)r2[1], (float)r2[2], (float)r2[3]));
}

static pxr::GfMatrix4d GetPrimLocalToParentTransform(const pxr::UsdPrim &prim,
                                                      pxr::UsdTimeCode timeCode = 0.0)
{
    pxr::GfMatrix4d localTransform(1.0);
    if (prim.IsA<pxr::UsdGeomXformable>())
    {
        pxr::UsdGeomXformable xformable(prim);
        bool resetsStack = false;
        xformable.GetLocalTransformation(&localTransform, &resetsStack, timeCode);
    }
    return localTransform;
}

static void TraversePrimToPrimLists(const pxr::UsdPrim &root, USDPrimLists *out)
{
    YBI_ASSERT(out);

    const pxr::Usd_PrimFlagsPredicate filterPredicate =
        pxr::UsdPrimIsActive && pxr::UsdPrimIsLoaded && !pxr::UsdPrimIsAbstract;

    std::vector<pxr::UsdPrim> stack;
    stack.push_back(root);

    while (!stack.empty())
    {
        const pxr::UsdPrim prim = stack.back();
        stack.pop_back();

        bool pushChildren = true;
        if (prim.IsInstance())
        {
            out->instances.push_back(prim);
            pushChildren = false;
        }
        else if (prim.IsA<pxr::UsdGeomPointInstancer>())
        {
            out->pointInstancers.push_back(prim);
            pushChildren = false;
        }
        else if (prim.IsA<pxr::UsdGeomMesh>())
        {
            out->meshes.push_back(prim);
        }
        else if (prim.IsA<pxr::UsdGeomBasisCurves>())
        {
            out->curves.push_back(prim);
        }

        if (!pushChildren)
        {
            continue;
        }

        std::vector<pxr::UsdPrim> children;
        for (const pxr::UsdPrim &child : prim.GetFilteredChildren(filterPredicate))
        {
            children.push_back(child);
        }
        for (size_t i = children.size(); i > 0; i--)
        {
            stack.push_back(children[i - 1]);
        }
    }
}

static void EnqueuePrototypePath(const std::string &prototypePath,
                                 std::unordered_map<std::string, int> *pathToPrototypeIndex,
                                 std::vector<std::string> *prototypePaths)
{
    YBI_ASSERT(pathToPrototypeIndex);
    YBI_ASSERT(prototypePaths);

    if (prototypePath.empty())
    {
        return;
    }
    if (pathToPrototypeIndex->find(prototypePath) != pathToPrototypeIndex->end())
    {
        return;
    }

    const int index = static_cast<int>(prototypePaths->size());
    pathToPrototypeIndex->emplace(prototypePath, index);
    prototypePaths->push_back(prototypePath);
}

static bool CollectPrototypeDependencies(const USDPrimLists &lists,
                                         std::unordered_map<std::string, int> *pathToPrototypeIndex,
                                         std::vector<std::string> *prototypePaths,
                                         std::string *error)
{
    for (const pxr::UsdPrim &instancePrim : lists.instances)
    {
        const pxr::UsdPrim prototype = instancePrim.GetPrototype();
        if (!prototype)
        {
            SetError(error, "instance without prototype at " + instancePrim.GetPath().GetString());
            return false;
        }
        EnqueuePrototypePath(prototype.GetPath().GetString(), pathToPrototypeIndex, prototypePaths);
    }

    for (const pxr::UsdPrim &pointInstancerPrim : lists.pointInstancers)
    {
        pxr::UsdGeomPointInstancer pointInstancer(pointInstancerPrim);
        pxr::SdfPathVector prototypePathsForInstancer;
        if (!pointInstancer.GetPrototypesRel().GetTargets(&prototypePathsForInstancer))
        {
            SetError(error,
                     "failed to read point instancer prototypes at " +
                         pointInstancerPrim.GetPath().GetString());
            return false;
        }

        for (const pxr::SdfPath &path : prototypePathsForInstancer)
        {
            EnqueuePrototypePath(path.GetString(), pathToPrototypeIndex, prototypePaths);
        }
    }

    return true;
}

static bool CollectPrototypePrimLists(const pxr::UsdStageRefPtr &stage,
                                      const USDPrimLists &rootPrimLists,
                                      std::unordered_map<std::string, int> *pathToPrototypeIndex,
                                      std::vector<std::string> *prototypePaths,
                                      std::vector<USDPrimLists> *prototypePrimLists,
                                      std::string *error)
{
    YBI_ASSERT(pathToPrototypeIndex);
    YBI_ASSERT(prototypePaths);
    YBI_ASSERT(prototypePrimLists);

    if (!CollectPrototypeDependencies(rootPrimLists, pathToPrototypeIndex, prototypePaths, error))
    {
        return false;
    }

    for (size_t prototypeIndex = 0; prototypeIndex < prototypePaths->size(); prototypeIndex++)
    {
        if (prototypePrimLists->size() < prototypePaths->size())
        {
            prototypePrimLists->resize(prototypePaths->size());
        }

        const std::string &prototypePath = (*prototypePaths)[prototypeIndex];
        const pxr::UsdPrim prototypePrim = stage->GetPrimAtPath(pxr::SdfPath(prototypePath));
        if (!prototypePrim)
        {
            SetError(error, "prototype prim not found: " + prototypePath);
            return false;
        }

        USDPrimLists primLists = {};
        primLists.ownerPath = prototypePath;
        TraversePrimToPrimLists(prototypePrim, &primLists);
        (*prototypePrimLists)[prototypeIndex] = std::move(primLists);

        if (!CollectPrototypeDependencies(
                (*prototypePrimLists)[prototypeIndex], pathToPrototypeIndex, prototypePaths, error))
        {
            return false;
        }
    }

    YBI_ASSERT(prototypePrimLists->size() == prototypePaths->size());
    return true;
}

static bool ResolvePrototypeChildSceneIndex(
    const std::string &prototypePath,
    const std::unordered_map<std::string, int> &pathToPrototypeIndex,
    const std::vector<std::unique_ptr<BuildScene>> &prototypeScenes,
    std::unordered_map<int, uint32_t> *childSceneIndexByPrototype,
    BuildScene *outScene,
    uint32_t *outChildSceneIndex,
    std::string *error)
{
    auto found = pathToPrototypeIndex.find(prototypePath);
    if (found == pathToPrototypeIndex.end())
    {
        SetError(error, "missing prototype path: " + prototypePath);
        return false;
    }

    const int prototypeIndex = found->second;
    if (prototypeIndex < 0 || prototypeIndex >= static_cast<int>(prototypeScenes.size()))
    {
        SetError(error, "prototype index out of range: " + prototypePath);
        return false;
    }
    if (!prototypeScenes[prototypeIndex])
    {
        SetError(error, "prototype scene missing when building parent: " + prototypePath);
        return false;
    }

    auto childFound = childSceneIndexByPrototype->find(prototypeIndex);
    if (childFound != childSceneIndexByPrototype->end())
    {
        *outChildSceneIndex = childFound->second;
        return true;
    }

    const uint32_t childSceneIndex = static_cast<uint32_t>(outScene->childScenes.size());
    outScene->childScenes.push_back(prototypeScenes[prototypeIndex].get());
    childSceneIndexByPrototype->emplace(prototypeIndex, childSceneIndex);
    *outChildSceneIndex = childSceneIndex;
    return true;
}

static bool AppendPointInstancerInstances(
    const pxr::UsdPrim &pointInstancerPrim,
    const std::unordered_map<std::string, int> &pathToPrototypeIndex,
    const std::vector<std::unique_ptr<BuildScene>> &prototypeScenes,
    std::unordered_map<int, uint32_t> *childSceneIndexByPrototype,
    pxr::UsdGeomXformCache *xformCache,
    BuildScene *outScene,
    std::string *error,
    pxr::UsdTimeCode timeCode = 0.0)
{
    pxr::UsdGeomPointInstancer pointInstancer(pointInstancerPrim);

    pxr::VtIntArray protoIndices;
    pxr::VtVec3fArray positions;
    pxr::VtQuatfArray orientations;
    pxr::VtVec3fArray scales;
    pxr::SdfPathVector prototypePaths;

    if (!pointInstancer.GetProtoIndicesAttr().Get(&protoIndices, timeCode) ||
        !pointInstancer.GetPositionsAttr().Get(&positions, timeCode) ||
        !pointInstancer.GetPrototypesRel().GetTargets(&prototypePaths))
    {
        SetError(error, "failed to read point instancer arrays at " +
                            pointInstancerPrim.GetPath().GetString());
        return false;
    }

    const pxr::UsdAttribute orientationsAttr = pointInstancer.GetOrientationsfAttr();
    const pxr::UsdAttribute scalesAttr = pointInstancer.GetScalesAttr();
    if (orientationsAttr.HasValue())
    {
        orientationsAttr.Get(&orientations, timeCode);
    }
    if (scalesAttr.HasValue())
    {
        scalesAttr.Get(&scales, timeCode);
    }

    if (protoIndices.size() != positions.size())
    {
        SetError(error, "point instancer protoIndices/positions size mismatch at " +
                            pointInstancerPrim.GetPath().GetString());
        return false;
    }
    if (!orientations.empty() && orientations.size() != positions.size())
    {
        SetError(error, "point instancer orientations/positions size mismatch at " +
                            pointInstancerPrim.GetPath().GetString());
        return false;
    }
    if (!scales.empty() && scales.size() != positions.size())
    {
        SetError(error, "point instancer scales/positions size mismatch at " +
                            pointInstancerPrim.GetPath().GetString());
        return false;
    }

    const pxr::GfMatrix4d pointInstancerLocalToWorld =
        xformCache->GetLocalToWorldTransform(pointInstancerPrim);
    pxr::UsdStageRefPtr stage = pointInstancerPrim.GetStage();

    for (size_t i = 0; i < protoIndices.size(); i++)
    {
        const int protoIndex = protoIndices[i];
        if (protoIndex < 0 || protoIndex >= static_cast<int>(prototypePaths.size()))
        {
            SetError(error, "point instancer prototype index out of range at " +
                                pointInstancerPrim.GetPath().GetString());
            return false;
        }

        const pxr::SdfPath prototypePath = prototypePaths[protoIndex];
        uint32_t childSceneIndex = 0;
        if (!ResolvePrototypeChildSceneIndex(prototypePath.GetString(),
                                             pathToPrototypeIndex,
                                             prototypeScenes,
                                             childSceneIndexByPrototype,
                                             outScene,
                                             &childSceneIndex,
                                             error))
        {
            return false;
        }

        const pxr::UsdPrim prototypePrim = stage->GetPrimAtPath(prototypePath);
        if (!prototypePrim)
        {
            SetError(error, "point instancer missing prototype prim: " + prototypePath.GetString());
            return false;
        }

        const pxr::GfMatrix4d prototypeLocalToParent =
            GetPrimLocalToParentTransform(prototypePrim, timeCode);

        pxr::GfMatrix4d scaleFromPrototype(1.0);
        if (!scales.empty())
        {
            scaleFromPrototype.SetScale(scales[i]);
        }

        pxr::GfMatrix4d orientationFromScale(1.0);
        if (!orientations.empty())
        {
            orientationFromScale.SetRotate(orientations[i]);
        }

        pxr::GfMatrix4d translationFromOrientation(1.0);
        translationFromOrientation.SetTranslateOnly(positions[i]);

        // Most-local (right) to least-local (left):
        // prototype local->parent, then S, R, T, then point instancer local->world.
        const pxr::GfMatrix4d worldFromLocal = pointInstancerLocalToWorld *
                                               translationFromOrientation * orientationFromScale *
                                               scaleFromPrototype * prototypeLocalToParent;

        SceneInstance instance = {};
        instance.worldFromLocal = ConvertAffineTransform(worldFromLocal);
        instance.childSceneIndex = childSceneIndex;
        outScene->instances.push_back(instance);
    }

    return true;
}

static bool BuildSceneFromPrimLists(const USDPrimLists &primLists,
                                    const std::unordered_map<std::string, int> &pathToPrototypeIndex,
                                    const std::vector<std::unique_ptr<BuildScene>> &prototypeScenes,
                                    pxr::UsdGeomXformCache *xformCache,
                                    BuildScene *outScene,
                                    std::string *error)
{
    YBI_ASSERT(xformCache);
    YBI_ASSERT(outScene);

    outScene->path = primLists.ownerPath;

    for (const pxr::UsdPrim &meshPrim : primLists.meshes)
    {
        outScene->meshes.push_back(
            {meshPrim.GetPath().GetString(),
             ConvertAffineTransform(xformCache->GetLocalToWorldTransform(meshPrim))});
    }

    for (const pxr::UsdPrim &curvePrim : primLists.curves)
    {
        outScene->curves.push_back(
            {curvePrim.GetPath().GetString(),
             ConvertAffineTransform(xformCache->GetLocalToWorldTransform(curvePrim))});
    }

    std::unordered_map<int, uint32_t> childSceneIndexByPrototype;

    for (const pxr::UsdPrim &instancePrim : primLists.instances)
    {
        const pxr::UsdPrim prototype = instancePrim.GetPrototype();
        if (!prototype)
        {
            SetError(error, "instance without prototype at " + instancePrim.GetPath().GetString());
            return false;
        }

        uint32_t childSceneIndex = 0;
        if (!ResolvePrototypeChildSceneIndex(prototype.GetPath().GetString(),
                                             pathToPrototypeIndex,
                                             prototypeScenes,
                                             &childSceneIndexByPrototype,
                                             outScene,
                                             &childSceneIndex,
                                             error))
        {
            return false;
        }

        SceneInstance instance = {};
        instance.worldFromLocal =
            ConvertAffineTransform(xformCache->GetLocalToWorldTransform(instancePrim));
        instance.childSceneIndex = childSceneIndex;
        outScene->instances.push_back(instance);
    }

    for (const pxr::UsdPrim &pointInstancerPrim : primLists.pointInstancers)
    {
        if (!AppendPointInstancerInstances(pointInstancerPrim,
                                           pathToPrototypeIndex,
                                           prototypeScenes,
                                           &childSceneIndexByPrototype,
                                           xformCache,
                                           outScene,
                                           error))
        {
            return false;
        }
    }

    return true;
}

} // namespace

bool BuildInstanceDAGFromUSD(const pxr::UsdStageRefPtr &stage, Scene *scene, std::string *error)
{
    if (!stage || !scene)
    {
        SetError(error, "invalid stage or scene");
        return false;
    }

    USDPrimLists rootPrimLists = {};
    rootPrimLists.ownerPath = stage->GetPseudoRoot().GetPath().GetString();
    TraversePrimToPrimLists(stage->GetPseudoRoot(), &rootPrimLists);

    std::unordered_map<std::string, int> pathToPrototypeIndex;
    std::vector<std::string> prototypePaths;
    std::vector<USDPrimLists> prototypePrimLists;
    if (!CollectPrototypePrimLists(stage,
                                   rootPrimLists,
                                   &pathToPrototypeIndex,
                                   &prototypePaths,
                                   &prototypePrimLists,
                                   error))
    {
        return false;
    }

    std::vector<int> reversedPrototypeOrder;
    reversedPrototypeOrder.reserve(prototypePaths.size());
    for (size_t i = 0; i < prototypePaths.size(); i++)
    {
        reversedPrototypeOrder.push_back(static_cast<int>(i));
    }
    std::reverse(reversedPrototypeOrder.begin(), reversedPrototypeOrder.end());

    BuildResult build = {};
    build.prototypeScenes.resize(prototypePrimLists.size());

    pxr::UsdGeomXformCache xformCache(0.0);

    for (int prototypeIndex : reversedPrototypeOrder)
    {
        std::unique_ptr<BuildScene> sceneForPrototype = std::make_unique<BuildScene>();
        if (!BuildSceneFromPrimLists(prototypePrimLists[(size_t)prototypeIndex],
                                     pathToPrototypeIndex,
                                     build.prototypeScenes,
                                     &xformCache,
                                     sceneForPrototype.get(),
                                     error))
        {
            return false;
        }
        build.prototypeScenes[(size_t)prototypeIndex] = std::move(sceneForPrototype);
    }

    build.rootScene = std::make_unique<BuildScene>();
    if (!BuildSceneFromPrimLists(rootPrimLists,
                                 pathToPrototypeIndex,
                                 build.prototypeScenes,
                                 &xformCache,
                                 build.rootScene.get(),
                                 error))
    {
        return false;
    }

    (void)build;
    return true;
}

YBI_NAMESPACE_END
