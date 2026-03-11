#pragma once

#include <string>
#include <vector>

namespace ybi
{

struct ScenePool;

struct USDLoadOptions
{
    std::vector<std::string> purposes = {"default", "render"};
};

void LoadUSDScene(ScenePool *scenePool, const std::string &filePath);
void LoadUSDScene(ScenePool *scenePool, const std::string &filePath, const USDLoadOptions &options);

} // namespace ybi
