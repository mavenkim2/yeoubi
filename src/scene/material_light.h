#pragma once

#include "util/float3.h"

#include <cstdint>
#include <string>
#include <vector>

YBI_NAMESPACE_BEGIN

enum MaterialFlags : uint32_t
{
    MATERIAL_FLAG_HAS_EMISSION = 1u << 0u,
};

struct PackedMaterial
{
    float3 baseColor = make_float3(0.7f, 0.7f, 0.7f);
    float roughness = 1.0f;
    float3 emissiveColor = make_float3(0.0f, 0.0f, 0.0f);
    float metallic = 0.0f;
    float ior = 1.5f;
    float opacity = 1.0f;
    uint32_t flags = 0u;
    uint32_t _padding0 = 0u;
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
    float3 position = make_float3(0.0f, 0.0f, 0.0f);
    float _padding0 = 0.0f;
    float3 direction = make_float3(0.0f, 0.0f, -1.0f);
    float _padding1 = 0.0f;
    float3 tangent = make_float3(1.0f, 0.0f, 0.0f);
    float _padding2 = 0.0f;
    float3 bitangent = make_float3(0.0f, 1.0f, 0.0f);
    float emissionScale = 0.0f;
    float3 color = make_float3(0.0f, 0.0f, 0.0f);
    float width = 0.0f;
    float height = 0.0f;
    float radius = 0.0f;
    float length = 0.0f;
    float areaScale = 0.0f;
    float _padding3 = 0.0f;
};

struct LightInfo
{
    std::string lightPath;
    std::string texturePath;
    PackedLight packed = {};
};

YBI_NAMESPACE_END
