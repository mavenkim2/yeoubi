#pragma once

#include <string>
#include <vector>

namespace ybi
{
namespace texture
{

bool LoadExrRgba(const std::string &path,
                 int *outWidth,
                 int *outHeight,
                 std::vector<float> *outRgba,
                 std::string *outReason,
                 bool flipVertical);

} // namespace texture
} // namespace ybi
