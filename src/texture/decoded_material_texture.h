#pragma once

#include "scene/scene.h"
#include "texture/texture_format.h"

#include <cstdint>
#include <string>
#include <vector>

namespace ybi
{

struct DecodedMaterialTexture
{
    bool valid = false;
    uint32_t udim = 1001u;
    int width = 0;
    int height = 0;
    TextureFormat format = TextureFormat::RGBA8_UNORM;
    std::vector<uint8_t> pixels;
    std::string sourcePath;
    std::string textureName;
    TextureWrapMode wrapS = TEXTURE_WRAP_MODE_REPEAT;
    TextureWrapMode wrapT = TEXTURE_WRAP_MODE_REPEAT;
};

} // namespace ybi
