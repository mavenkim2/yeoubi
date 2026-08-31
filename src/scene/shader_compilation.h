#include "render/virtual_texture_manager.h"
#include "scene/material.h"

#include <unordered_map>

namespace ybi
{
int CompileShader(
    const MaterialGraph &graph,
    MaterialNodeId outputNode,
    const std::unordered_map<const MaterialNode *, MaterialTextureHandle> &textureNodeToHandle,
    std::vector<uint32_t> &nodes);
} // namespace ybi
