#pragma once

#include <string>

YBI_NAMESPACE_BEGIN

struct ScenePool;
void LoadUSDScene(ScenePool *scenePool, const std::string &filePath);

YBI_NAMESPACE_END
