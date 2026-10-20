#pragma once

#include "core/affine_transform.h"
#include "core/vec3f.h"
#include <cstdint>
#include <string>
#include <vector>

namespace ybi
{

struct EnvmapLight
{
    std::string texturePath;
    std::string textureFormat;
    AffineTransform localToWorld;
    Vec3f radianceScale;

    uint32_t width = 0;
    uint32_t height = 0;
    std::vector<Vec3f> pixels; // Linear RGB radiance, including radianceScale.
    std::vector<float> conditionalCDF;
    std::vector<float> marginalCDF;
    float totalWeight = 0.f;
};

bool CreateEnvmapCDF(EnvmapLight &light);

} // namespace ybi
