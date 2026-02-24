#pragma once

#include <string>
#include <vector>

namespace ybi
{
namespace usd_ntc
{

bool LoadExrRgba(const std::string &path, int &width, int &height, std::vector<float> &rgba, std::string &reason);

} // namespace usd_ntc
} // namespace ybi
