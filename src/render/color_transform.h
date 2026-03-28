#pragma once

#include <cstdint>
#include <string>
#include <vector>

namespace ybi
{
namespace render
{

bool ApplyAcesSdrDisplayTransform(const std::vector<float> &linearRgb,
                                  int width,
                                  int height,
                                  std::vector<uint8_t> *outRgba,
                                  std::string *outError);

} // namespace render
} // namespace ybi
