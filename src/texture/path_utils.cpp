#include "path_utils.h"

#include <cctype>
#include <filesystem>

namespace ybi
{
namespace texture
{

std::string LowerExt(const std::string &path)
{
    std::string ext = std::filesystem::path(path).extension().string();
    for (char &c : ext)
    {
        c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    }
    return ext;
}

} // namespace texture
} // namespace ybi
