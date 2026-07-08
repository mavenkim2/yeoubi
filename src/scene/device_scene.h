#pragma once

#include "core/vec2f.h"
#include "core/vec3f.h"
#include "device/device_memory.h"

namespace ybi
{

enum class AttributeType
{
    None,
    Vertex,
    FaceVarying,
};

struct DeviceMaterial
{
    // NOTE: = ~0u if no udims
    uint32_t udimBase;
    uint32_t udimTableOffset;
};

struct Object
{
    AttributeType uvType;
    uint32_t uvOffset;
    uint32_t positionOffset;
    uint32_t triangleIndexOffset;
    uint32_t triangleMaterialOffset;

    int materialIndex;
};

struct DeviceScene
{
    Object *objects;
    DeviceMaterial *materials;

    Vec3f *positions;
    Vec2f *float2s;
    int *triangleIndices;
    int *triangleMaterialIndices;
};

} // namespace ybi
