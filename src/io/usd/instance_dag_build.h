#pragma once

#include <pxr/usd/usd/prim.h>
#include <pxr/usd/usd/stage.h>
#include "util/vec3.h"
#include "util/float3x4.h"
#include <string>
#include <unordered_map>
#include <vector>

namespace ybi
{

struct USDPrimLists
{
    std::string ownerPath;
    std::vector<pxr::UsdPrim> meshes;
    std::vector<pxr::UsdPrim> curves;
    std::vector<pxr::UsdPrim> instances;
    std::vector<pxr::UsdPrim> pointInstancers;
};

struct USDPrototypePrimLists
{
    std::unordered_map<std::string, int> prototypePathToIndex;
    std::vector<std::string> prototypePaths;
    std::vector<USDPrimLists> prototypePrimLists;
};

struct USDBuildSceneMesh
{
    std::string path;
    Float3x4 parentFromLocal;
};

struct USDBuildSceneCurve
{
    std::string path;
    Float3x4 parentFromLocal;
};

struct USDBuildSceneInstance
{
    Float3x4 parentFromLocal;
    uint32_t childSceneIndex;
};

struct USDBuildScene
{
    std::string path;
    std::vector<USDBuildSceneMesh> meshes;
    std::vector<USDBuildSceneCurve> curves;
    std::vector<USDBuildSceneInstance> instances;
};

struct USDBuildSceneDAG
{
    std::vector<USDBuildScene> scenes;
    uint32_t rootSceneIndex;
};

void TraverseUSDPrimLists(const pxr::UsdPrim &root, USDPrimLists *out);

bool CollectUSDPrototypePrimLists(const pxr::UsdStageRefPtr &stage,
                                  const USDPrimLists &rootPrimLists,
                                  USDPrototypePrimLists *out,
                                  std::string *error);

bool BuildInstanceDAGFromUSD(const pxr::UsdStageRefPtr &stage,
                             USDBuildSceneDAG *out,
                             std::string *error);

} // namespace ybi
