#pragma once

#include <pxr/usd/usd/stage.h>
#include <string>

YBI_NAMESPACE_BEGIN

struct Scene;

// Builds scene->instances / scene->sceneRefs / scene->rootRefs from USD instances and point
// instancers. Prototypes are resolved to SceneRefRange (not Scene instances).
bool BuildInstanceDAGFromUSD(const pxr::UsdStageRefPtr &stage, Scene *scene, std::string *error);

YBI_NAMESPACE_END
