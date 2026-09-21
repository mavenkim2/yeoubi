#pragma once

#include "render/virtual_texture_manager.h"
#include "scene/device_scene.h"
#include "scene/material.h"

#include <unordered_map>

namespace ybi
{

DeviceMaterialExpression CompileShader(
    const MaterialGraph &graph,
    MaterialNodeId outputNode,
    const std::unordered_map<const MaterialNode *, MaterialTextureHandle> &textureNodeToHandle,
    const std::unordered_map<std::string, AttributeId> &attributeIdMap,
    std::vector<uint32_t> &nodes);

} // namespace ybi
