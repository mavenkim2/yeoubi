#pragma once

#include "scene/scene.h"

#include <cstdint>

namespace ybi
{

uint64_t ComputeMeshHostBytes(const Mesh &mesh);
uint64_t ComputeMeshVectorHostBytes(const std::vector<Mesh> &meshes);
uint64_t ComputeSceneHostBytes(const Scene &scene);
uint64_t ComputeScenePoolHostBytes(const ScenePool &scenePool);

} // namespace ybi
