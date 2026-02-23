#pragma once

#include "scene/scene.h"

#include <string>
#include <vector>

namespace testbvh
{

struct DecodedMaterialTexture
{
    bool valid = false;
    int width = 0;
    int height = 0;
    std::vector<unsigned char> rgba8;
    std::string ntcPath;
    std::string textureName;
};

bool DecodeNtcDiffuseTextures(const std::vector<ybi::ScenePool::MaterialInfo> &materials,
                              const std::string &ntcDir,
                              std::vector<DecodedMaterialTexture> *outTextures,
                              std::string *outError);

} // namespace testbvh
