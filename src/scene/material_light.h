#pragma once

#include "util/float3x4.h"
#include "util/vec3.h"

#include <cstdint>
#include <string>
#include <vector>

namespace ybi
{

enum MaterialFlags : uint32_t
{
    MATERIAL_FLAG_HAS_EMISSION = 1u << 0u,
};

struct PackedMaterial
{
    Vec3 baseColor = Vec3(0.18f, 0.18f, 0.18f);
    float roughness = 0.5f;
    Vec3 emissiveColor = Vec3(0.0f, 0.0f, 0.0f);
    float metallic = 0.0f;
    Vec3 specularColor = Vec3(0.0f, 0.0f, 0.0f);
    float ior = 1.5f;
    float opacity = 1.0f;
    float clearcoat = 0.0f;
    float clearcoatRoughness = 0.01f;
    float opacityThreshold = 0.0f;
    uint32_t flags = 0u;
    uint32_t useSpecularWorkflow = 0u;
    uint32_t _padding0 = 0u;
    uint32_t _padding1 = 0u;
};

enum class LightType : uint32_t
{
    Dome = 0,
    Distant = 1,
    Rect = 2,
    Disk = 3,
    Sphere = 4,
    Cylinder = 5,
};

enum LightFlags : uint32_t
{
    LIGHT_FLAG_NORMALIZED = 1u << 0u,
    LIGHT_FLAG_ONE_SIDED = 1u << 1u,
};

struct PackedLight
{
    uint32_t type = 0u;
    uint32_t flags = 0u;
    Float3x4 worldFromLocal = Float3x4::Identity();
    float emissionScale = 0.0f;
    Vec3 color = Vec3(0.0f, 0.0f, 0.0f);
    float angleRadians = 0.0f;
    float width = 0.0f;
    float height = 0.0f;
    float radius = 0.0f;
    float length = 0.0f;
    float cosThetaMax = 1.0f;
    float areaScale = 0.0f;
    float solidAngle = 0.0f;
    float selectionWeight = 1.0f;
    uint32_t shadowExcludeOffset = 0u;
    uint32_t shadowExcludeCount = 0u;
};

struct LightInfo
{
    std::string lightPath;
    std::string texturePath;
    std::vector<std::string> shadowExcludePaths;
    PackedLight packed = {};
};

} // namespace ybi
