#pragma once

#include "render/virtual_texture_interop.h"
#include "scene/device_scene.h"
#include "scene/material.h"

#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

namespace ybi
{

struct Device;

struct MaterialCompilationResult
{
    std::string source;
    std::vector<DeviceMaterial3> materials;
    bool succeeded = false;
};

class MaterialCompiler
{
    uint32_t numMaterials = 0;
    std::vector<std::string> generatedMaterialSources;

public:
    std::string Compile(const MaterialGraph &graph);
    std::string Compile(const MaterialGraph &graph,
                        const std::unordered_map<const MaterialNode *, MaterialTextureHandle>
                            &textureNodeToHandle);
    MaterialCompilationResult Finalize(Device &device) const;
};

} // namespace ybi
