#pragma once

#include "scene/scene.h"

#include <cuda_runtime_api.h>

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

struct MaterialTextureRef
{
    unsigned long long textureObject = 0ull;
    int width = 0;
    int height = 0;
    int valid = 0;
};

struct UploadedMaterialTextures
{
    std::vector<MaterialTextureRef> refs;
    std::vector<cudaArray_t> arrays;
    std::vector<cudaTextureObject_t> textureObjects;
};

bool DecodeNtcDiffuseTextures(const std::vector<ybi::ScenePool::MaterialInfo> &materials,
                              const std::string &ntcDir,
                              std::vector<DecodedMaterialTexture> *outTextures,
                              std::string *outError);

bool UploadDecodedTexturesToCuda(const std::vector<DecodedMaterialTexture> &decodedTextures,
                                 UploadedMaterialTextures *outTextures,
                                 std::string *outError);

void DestroyUploadedTextures(UploadedMaterialTextures *textures);

} // namespace testbvh
