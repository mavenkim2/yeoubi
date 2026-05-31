#pragma once

#include <cstdint>
#include <string>
#include <utility>
#include <vector>

namespace ybi
{
namespace texture
{

bool TryFindUdimDigits(const std::string &path, uint32_t &udim, size_t &digitPos);
std::string StripUdimFromPath(const std::string &path);
bool CollectUdimPaths(const std::string &path,
                      std::vector<std::pair<uint32_t, std::string>> &out,
                      std::string &reason);

} // namespace texture
} // namespace ybi
