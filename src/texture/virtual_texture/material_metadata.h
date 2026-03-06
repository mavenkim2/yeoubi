#pragma once

#include "render/launch_params.h"
#include "scene/scene.h"
#include "texture/virtual_texture/manager.h"

#include <cstddef>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

namespace ybi
{
namespace texture
{

struct VirtualTextureMaterialSource
{
    uint32_t materialIndex = 0u;
    std::string texturePath;
    TextureWrapMode wrapS = TEXTURE_WRAP_MODE_UNKNOWN;
    TextureWrapMode wrapT = TEXTURE_WRAP_MODE_UNKNOWN;
};

struct VirtualTextureMaterialBuildResult
{
    std::vector<LaunchParams::MaterialTextureRef> materialTextureRefs;
    std::vector<VirtualTextureRegisterInput> registrations;
    uint64_t totalVirtualTextureBytes = 0u;
    uint32_t mappedMaterialCount = 0u;
    uint32_t activeUdimCount = 0u;
};

bool BuildVirtualTextureMaterialMetadata(
    size_t materialCount,
    int semanticIndex,
    int materialTextureRefSemanticCount,
    int materialTextureRefStride,
    const std::vector<VirtualTextureMaterialSource> &sources,
    const std::unordered_map<unsigned int, std::string> &tileFiles,
    VirtualTextureMaterialBuildResult *out,
    std::string *outError);

} // namespace texture
} // namespace ybi
