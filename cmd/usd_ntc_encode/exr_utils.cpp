#include "exr_utils.h"
#include "texture/exr_io.h"

#include <string>
#include <vector>

namespace ybi
{
namespace usd_ntc
{

bool LoadExrRgba(const std::string &path, int &width, int &height, std::vector<float> &rgba, std::string &reason)
{
    return ybi::texture::LoadExrRgba(path, &width, &height, &rgba, &reason, false);
}

} // namespace usd_ntc
} // namespace ybi
