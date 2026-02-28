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

struct SceneInstance
{
    float3x4 parentFromLocal;
    uint32_t childSceneIndex;
};

struct SceneMesh
{
    std::string path;
    float3x4 parentFromLocal;
};

struct SceneCurve
{
    std::string path;
    float3x4 parentFromLocal;
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

static void
TraversePrimToPrimListsImpl(const pxr::UsdPrim &root, USDPrimLists *out, bool skipPrototypePrims)
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

        if (skipPrototypePrims && (prim.IsInPrototype() || prim.IsInstanceProxy()))
        {
            continue;
        }

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

static bool
CollectPrototypeDependencies(const USDPrimLists &lists,
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
        EnqueuePrototypePath(
            prototype.GetPath().GetString(), pathToPrototypeIndex, prototypePaths);
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

static bool
CollectPrototypePrimListsImpl(const pxr::UsdStageRefPtr &stage,
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
        TraversePrimToPrimListsImpl(prototypePrim, &primLists, /*skipPrototypePrims*/ false);
        (*prototypePrimLists)[prototypeIndex] = std::move(primLists);

        if (!CollectPrototypeDependencies((*prototypePrimLists)[prototypeIndex],
                                          pathToPrototypeIndex,
                                          prototypePaths,
                                          error))
        {
            return false;
        }
    }

    YBI_ASSERT(prototypePrimLists->size() == prototypePaths->size());
    return true;
}

static bool
GetPrototypeTopologicalOrder(const std::unordered_map<std::string, int> &pathToPrototypeIndex,
                             const std::vector<USDPrimLists> &prototypePrimLists,
                             std::vector<int> *outSortedOrder,
                             std::string *error)
{
    YBI_ASSERT(outSortedOrder);

    const int numPrototypes = static_cast<int>(prototypePrimLists.size());
    outSortedOrder->clear();
    outSortedOrder->reserve((size_t)numPrototypes);
    if (numPrototypes == 0)
    {
        return true;
    }

    // TODO: jagged array instead?
    std::vector<std::vector<int>> dependentsByPrototype((size_t)numPrototypes);
    std::vector<int> dependencyCount((size_t)numPrototypes, 0);

    for (int prototypeIndex = 0; prototypeIndex < numPrototypes; prototypeIndex++)
    {
        std::vector<uint8_t> hasDependency((size_t)numPrototypes, 0);
        const USDPrimLists &lists = prototypePrimLists[(size_t)prototypeIndex];

        for (const pxr::UsdPrim &instancePrim : lists.instances)
        {
            const pxr::UsdPrim dependencyPrototype = instancePrim.GetPrototype();
            if (!dependencyPrototype)
            {
                SetError(error,
                         "instance without prototype at " + instancePrim.GetPath().GetString());
                return false;
            }

            auto depFound = pathToPrototypeIndex.find(dependencyPrototype.GetPath().GetString());
            if (depFound == pathToPrototypeIndex.end())
            {
                SetError(error,
                         "prototype dependency missing from index map: " +
                             dependencyPrototype.GetPath().GetString());
                return false;
            }

            const int dependencyIndex = depFound->second;
            if (dependencyIndex < 0 || dependencyIndex >= numPrototypes)
            {
                SetError(error,
                         "prototype dependency index out of range: " +
                             dependencyPrototype.GetPath().GetString());
                return false;
            }
            if (!hasDependency[(size_t)dependencyIndex])
            {
                hasDependency[(size_t)dependencyIndex] = 1;
                dependencyCount[(size_t)prototypeIndex]++;
                dependentsByPrototype[(size_t)dependencyIndex].push_back(prototypeIndex);
            }
        }

        for (const pxr::UsdPrim &pointInstancerPrim : lists.pointInstancers)
        {
            pxr::UsdGeomPointInstancer pointInstancer(pointInstancerPrim);
            pxr::SdfPathVector dependencyPrototypePaths;
            if (!pointInstancer.GetPrototypesRel().GetTargets(&dependencyPrototypePaths))
            {
                SetError(error,
                         "failed to read point instancer prototypes at " +
                             pointInstancerPrim.GetPath().GetString());
                return false;
            }

            for (const pxr::SdfPath &dependencyPath : dependencyPrototypePaths)
            {
                auto depFound = pathToPrototypeIndex.find(dependencyPath.GetString());
                if (depFound == pathToPrototypeIndex.end())
                {
                    SetError(error,
                             "point instancer prototype dependency missing from index map: " +
                                 dependencyPath.GetString());
                    return false;
                }

                const int dependencyIndex = depFound->second;
                if (dependencyIndex < 0 || dependencyIndex >= numPrototypes)
                {
                    SetError(error,
                             "point instancer dependency index out of range: " +
                                 dependencyPath.GetString());
                    return false;
                }
                if (!hasDependency[(size_t)dependencyIndex])
                {
                    hasDependency[(size_t)dependencyIndex] = 1;
                    dependencyCount[(size_t)prototypeIndex]++;
                    dependentsByPrototype[(size_t)dependencyIndex].push_back(prototypeIndex);
                }
            }
        }
    }

    std::vector<int> ready;
    ready.reserve((size_t)numPrototypes);
    for (int prototypeIndex = 0; prototypeIndex < numPrototypes; prototypeIndex++)
    {
        if (dependencyCount[(size_t)prototypeIndex] == 0)
        {
            ready.push_back(prototypeIndex);
        }
    }

    std::vector<int> sortedOrder;
    sortedOrder.reserve((size_t)numPrototypes);

    size_t readyCursor = 0;
    while (readyCursor < ready.size())
    {
        const int next = ready[readyCursor++];
        outSortedOrder->push_back(next);

        for (int dependentIndex : dependentsByPrototype[(size_t)next])
        {
            int &count = dependencyCount[(size_t)dependentIndex];
            YBI_ASSERT(count > 0);
            count--;
            if (count == 0)
            {
                ready.push_back(dependentIndex);
            }
        }
    }

    if (outSortedOrder->size() != (size_t)numPrototypes)
    {
        SetError(error, "prototype dependency cycle detected while topologically sorting");
        return false;
    }

    return true;
}

static bool
ResolvePrototypeChildSceneIndex(const std::string &prototypePath,
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

static bool
AppendPointInstancerInstances(const pxr::UsdPrim &pointInstancerPrim,
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
        SetError(error,
                 "failed to read point instancer arrays at " +
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
        SetError(error,
                 "point instancer protoIndices/positions size mismatch at " +
                     pointInstancerPrim.GetPath().GetString());
        return false;
    }
    if (!orientations.empty() && orientations.size() != positions.size())
    {
        SetError(error,
                 "point instancer orientations/positions size mismatch at " +
                     pointInstancerPrim.GetPath().GetString());
        return false;
    }
    if (!scales.empty() && scales.size() != positions.size())
    {
        SetError(error,
                 "point instancer scales/positions size mismatch at " +
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
            SetError(error,
                     "point instancer prototype index out of range at " +
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
            SetError(error,
                     "point instancer missing prototype prim: " + prototypePath.GetString());
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
        const pxr::GfMatrix4d parentFromLocal = pointInstancerLocalToWorld *
                                                translationFromOrientation * orientationFromScale *
                                                scaleFromPrototype * prototypeLocalToParent;

        SceneInstance instance = {};
        instance.parentFromLocal = ConvertAffineTransform(parentFromLocal);
        instance.childSceneIndex = childSceneIndex;
        outScene->instances.push_back(instance);
    }

    return true;
}

static bool
BuildSceneFromPrimLists(const USDPrimLists &primLists,
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
        instance.parentFromLocal =
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

static bool
ExportBuildSceneDAG(const BuildResult &build, USDBuildSceneDAG *out, std::string *error)
{
    YBI_ASSERT(out);
    if (!build.rootScene)
    {
        SetError(error, "missing root build scene");
        return false;
    }

    const uint32_t numPrototypeScenes = static_cast<uint32_t>(build.prototypeScenes.size());
    const uint32_t rootSceneIndex = numPrototypeScenes;
    out->scenes.clear();
    out->scenes.resize(numPrototypeScenes + 1);
    out->rootSceneIndex = rootSceneIndex;

    std::unordered_map<const BuildScene *, uint32_t> sceneIndexByPtr;
    for (uint32_t i = 0; i < numPrototypeScenes; i++)
    {
        if (!build.prototypeScenes[i])
        {
            SetError(error, "missing prototype build scene");
            return false;
        }
        sceneIndexByPtr.emplace(build.prototypeScenes[i].get(), i);
    }
    sceneIndexByPtr.emplace(build.rootScene.get(), rootSceneIndex);

    auto copyScene = [&](const BuildScene *src, USDBuildScene *dst) -> bool {
        dst->path = src->path;
        dst->meshes.clear();
        dst->curves.clear();
        dst->instances.clear();

        dst->meshes.reserve(src->meshes.size());
        for (const SceneMesh &mesh : src->meshes)
        {
            dst->meshes.push_back({mesh.path, mesh.parentFromLocal});
        }

        dst->curves.reserve(src->curves.size());
        for (const SceneCurve &curve : src->curves)
        {
            dst->curves.push_back({curve.path, curve.parentFromLocal});
        }

        dst->instances.reserve(src->instances.size());
        for (const SceneInstance &instance : src->instances)
        {
            if (instance.childSceneIndex >= src->childScenes.size())
            {
                SetError(error, "instance child index out of range while exporting build scene");
                return false;
            }

            const BuildScene *childScene = src->childScenes[instance.childSceneIndex];
            auto childIndexIter = sceneIndexByPtr.find(childScene);
            if (childIndexIter == sceneIndexByPtr.end())
            {
                SetError(error, "child scene pointer missing while exporting build scene");
                return false;
            }

            dst->instances.push_back({instance.parentFromLocal, childIndexIter->second});
        }

        return true;
    };

    for (uint32_t i = 0; i < numPrototypeScenes; i++)
    {
        if (!copyScene(build.prototypeScenes[i].get(), &out->scenes[i]))
        {
            return false;
        }
    }

    if (!copyScene(build.rootScene.get(), &out->scenes[rootSceneIndex]))
    {
        return false;
    }

    return true;
}

} // namespace

void TraverseUSDPrimLists(const pxr::UsdPrim &root, USDPrimLists *out)
{
    if (!out)
    {
        return;
    }

    out->ownerPath = root.GetPath().GetString();
    out->meshes.clear();
    out->curves.clear();
    out->instances.clear();
    out->pointInstancers.clear();
    TraversePrimToPrimListsImpl(root, out, /*skipPrototypePrims*/ true);
}

bool CollectUSDPrototypePrimLists(const pxr::UsdStageRefPtr &stage,
                                  const USDPrimLists &rootPrimLists,
                                  USDPrototypePrimLists *out,
                                  std::string *error)
{
    if (!stage || !out)
    {
        SetError(error, "invalid stage or output");
        return false;
    }

    out->prototypePathToIndex.clear();
    out->prototypePaths.clear();
    out->prototypePrimLists.clear();

    return CollectPrototypePrimListsImpl(stage,
                                         rootPrimLists,
                                         &out->prototypePathToIndex,
                                         &out->prototypePaths,
                                         &out->prototypePrimLists,
                                         error);
}

bool BuildInstanceDAGFromUSD(const pxr::UsdStageRefPtr &stage,
                             USDBuildSceneDAG *out,
                             std::string *error)
{
    if (!stage || !out)
    {
        SetError(error, "invalid stage or output");
        return false;
    }

    USDPrimLists rootPrimLists = {};
    TraverseUSDPrimLists(stage->GetPseudoRoot(), &rootPrimLists);

    USDPrototypePrimLists prototypeData = {};
    if (!CollectUSDPrototypePrimLists(stage, rootPrimLists, &prototypeData, error))
    {
        return false;
    }

    BuildResult build = {};
    build.prototypeScenes.resize(prototypeData.prototypePrimLists.size());

    std::vector<int> prototypeBuildOrder;
    if (!GetPrototypeTopologicalOrder(prototypeData.prototypePathToIndex,
                                      prototypeData.prototypePrimLists,
                                      &prototypeBuildOrder,
                                      error))
    {
        return false;
    }

    pxr::UsdGeomXformCache xformCache(0.0);

    for (int prototypeIndex : prototypeBuildOrder)
    {
        std::unique_ptr<BuildScene> sceneForPrototype = std::make_unique<BuildScene>();
        if (!BuildSceneFromPrimLists(prototypeData.prototypePrimLists[(size_t)prototypeIndex],
                                     prototypeData.prototypePathToIndex,
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
                                 prototypeData.prototypePathToIndex,
                                 build.prototypeScenes,
                                 &xformCache,
                                 build.rootScene.get(),
                                 error))
    {
        return false;
    }

    return ExportBuildSceneDAG(build, out, error);
}

YBI_NAMESPACE_END
