#pragma once

#include "scene/scene.h"
#include "texture/decoded_material_texture.h"

#include <cuda_runtime_api.h>

#include <string>
#include <vector>

namespace testbvh
{

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

bool DecodeNtcDiffuseTextures(const std::vector<ybi::MaterialInfo> &materials,
                              std::vector<ybi::DecodedMaterialTexture> *outTextures,
                              std::string *outError);

bool UploadDecodedTexturesToCuda(const std::vector<ybi::DecodedMaterialTexture> &decodedTextures,
                                 UploadedMaterialTextures *outTextures,
                                 std::string *outError);

void DestroyUploadedTextures(UploadedMaterialTextures *textures);

} // namespace testbvh
