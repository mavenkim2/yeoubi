#pragma once

#include "texture/decoded_material_texture.h"

#include <string>
#include <vector>

namespace ybi
{
namespace lights
{

bool BuildDomeTextureCdf(const DecodedMaterialTexture &domeTexture,
                         std::vector<float> *outConditional,
                         std::vector<float> *outMarginal,
                         std::string *outError);

} // namespace lights
} // namespace ybi
