#pragma once

#include "bvh/usd_camera_utils.h"
#include "scene/subdivision_mesh.h"

#include <filesystem>
#include <string>

namespace ybi::testio
{
bool LoadSelectedSubdivFromJson(const std::filesystem::path &jsonPath,
                                ybi::SubdivisionMesh &meshOut,
                                UsdCameraInfo &usdCameraInfoOut,
                                std::string *outSubdivisionScheme = nullptr,
                                std::string *outCreasingMethod = nullptr,
                                std::string *outTriangleSubdivision = nullptr);
} // namespace ybi::testio
