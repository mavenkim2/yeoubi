#pragma once

#include <string>
#include <vector>

YBI_NAMESPACE_BEGIN

struct ScenePool;

struct USDLoadOptions
{
    std::vector<std::string> purposes = {"default", "render"};
};

void LoadUSDScene(ScenePool *scenePool, const std::string &filePath);
void LoadUSDScene(ScenePool *scenePool, const std::string &filePath, const USDLoadOptions &options);

YBI_NAMESPACE_END
