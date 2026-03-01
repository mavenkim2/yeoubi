#pragma once

#include <cstdint>
#include <string>
#include <unordered_map>

namespace ybi
{
namespace usd_ntc
{

bool TryFindUdimDigits(const std::string &path, uint32_t &udim, size_t &digitPos);
std::string StripUdimFromPath(const std::string &path);
bool CollectUdimPaths(const std::string &path, std::unordered_map<uint32_t, std::string> &out, std::string &reason);

} // namespace usd_ntc
} // namespace ybi
