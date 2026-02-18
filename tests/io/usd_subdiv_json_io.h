#pragma once

#include "bvh/usd_camera_utils.h"
#include "bvh/usd_subdiv_select.h"

#include <filesystem>

namespace ybi::testio
{
bool LoadSelectedSubdivFromJson(const std::filesystem::path &jsonPath,
                                SelectedSubdivMesh &selectedOut,
                                UsdCameraInfo &usdCameraInfoOut);
} // namespace ybi::testio
