#pragma once

#include "core/affine_transform.h"
#include "core/vec3f.h"
#include "core/vec4f.h"
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
    std::vector<Vec4f> pixels; // Linear RGBA; RGB includes radianceScale.
    std::vector<float> conditionalCDF;
    std::vector<float> marginalCDF;
    float totalWeight = 0.f;
};

bool LoadEnvmap(EnvmapLight &light, std::string *error);
void PrepareEnvmap(EnvmapLight &light);

} // namespace ybi
