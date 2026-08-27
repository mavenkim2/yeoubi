#include "scene/material.h"
#include "render/virtual_texture_manager.h"

#include <unordered_map>

namespace ybi
{
void CompileShader(const MaterialGraph &graph,
                   MaterialNodeId outputNode,
                   const std::unordered_map<const MaterialNode *, MaterialTextureHandle>
                       &textureNodeToHandle,
                   std::vector<uint32_t> &nodes);
} // namespace ybi
