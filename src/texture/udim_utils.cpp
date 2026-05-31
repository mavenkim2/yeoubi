#include "udim_utils.h"

#include <cctype>
#include <cstdlib>
#include <filesystem>
#include <string>
#include <utility>
#include <vector>

namespace ybi
{
namespace texture
{

static constexpr uint32_t kUdimMin = 1001;
static constexpr uint32_t kUdimMax = 1100;

namespace
{

bool FindUdimToken(const std::string &path, size_t &pos, size_t &len)
{
    pos = path.find("<UDIM>");
    len = 6;
    if (pos != std::string::npos)
    {
        return true;
    }
    pos = path.find("<udim>");
    len = 6;
    return pos != std::string::npos;
}

} // namespace

bool TryFindUdimDigits(const std::string &path, uint32_t &udim, size_t &digitPos)
{
    const size_t ext = path.rfind('.');
    if (ext == std::string::npos || ext < 4)
    {
        return false;
    }
    digitPos = ext - 4;
    for (size_t i = 0; i < 4; ++i)
    {
        if (!std::isdigit(static_cast<unsigned char>(path[digitPos + i])))
        {
            return false;
        }
    }
    if (digitPos > 0)
    {
        const char c = path[digitPos - 1];
        if (!(c == '.' || c == '_' || c == '-'))
        {
            return false;
        }
    }
    udim = static_cast<uint32_t>(std::strtoul(path.substr(digitPos, 4).c_str(), nullptr, 10));
    return udim >= kUdimMin && udim <= kUdimMax;
}

std::string StripUdimFromPath(const std::string &path)
{
    size_t tokenPos = std::string::npos;
    size_t tokenLen = 0;
    if (FindUdimToken(path, tokenPos, tokenLen))
    {
        std::string out = path;
        out.erase(tokenPos, tokenLen);
        if (tokenPos > 0 && tokenPos < out.size() && out[tokenPos - 1] == '.' && out[tokenPos] == '.')
        {
            out.erase(tokenPos - 1, 1);
        }
        if (tokenPos > 0 && tokenPos < out.size() &&
            (out[tokenPos - 1] == '_' || out[tokenPos - 1] == '-') && out[tokenPos] == '.')
        {
            out.erase(tokenPos - 1, 1);
        }
        return out;
    }

    uint32_t udim = 0;
    size_t digitPos = 0;
    if (TryFindUdimDigits(path, udim, digitPos))
    {
        std::string out = path;
        out.erase(digitPos, 4);
        if (digitPos > 0 && digitPos < out.size() && out[digitPos - 1] == '.' && out[digitPos] == '.')
        {
            out.erase(digitPos - 1, 1);
        }
        if (digitPos > 0 && digitPos < out.size() &&
            (out[digitPos - 1] == '_' || out[digitPos - 1] == '-') && out[digitPos] == '.')
        {
            out.erase(digitPos - 1, 1);
        }
        return out;
    }

    return path;
}

bool CollectUdimPaths(const std::string &path,
                      std::vector<std::pair<uint32_t, std::string>> &out,
                      std::string &reason)
{
    out.clear();

    size_t tokenPos = std::string::npos;
    size_t tokenLen = 0;
    if (FindUdimToken(path, tokenPos, tokenLen))
    {
        for (uint32_t udim = kUdimMin; udim <= kUdimMax; ++udim)
        {
            std::string candidate = path;
            candidate.replace(tokenPos, tokenLen, std::to_string(udim));
            if (std::filesystem::exists(candidate))
            {
                out.emplace_back(udim, std::move(candidate));
            }
        }
        if (out.empty())
        {
            reason = "no UDIM files found for " + path;
            return false;
        }
        return true;
    }

    uint32_t explicitUdim = 0;
    size_t digitPos = 0;
    if (TryFindUdimDigits(path, explicitUdim, digitPos))
    {
        const std::string prefix = path.substr(0, digitPos);
        const std::string suffix = path.substr(digitPos + 4);
        for (uint32_t udim = kUdimMin; udim <= kUdimMax; ++udim)
        {
            const std::string candidate = prefix + std::to_string(udim) + suffix;
            if (std::filesystem::exists(candidate))
            {
                out.emplace_back(udim, candidate);
            }
        }
        if (!out.empty())
        {
            return true;
        }
        if (std::filesystem::exists(path))
        {
            out.emplace_back(explicitUdim, path);
            return true;
        }
        reason = "missing texture file " + path;
        return false;
    }

    if (std::filesystem::exists(path))
    {
        out.emplace_back(kUdimMin, path);
        return true;
    }

    reason = "missing texture file " + path;
    return false;
}

} // namespace texture
} // namespace ybi
