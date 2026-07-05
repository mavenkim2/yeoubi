#pragma once

#include "device/device_memory.h"

namespace ybi
{

enum class AttributeType
{
    None,
    Vertex,
    FaceVarying,
};

struct Object
{
    AttributeType uvType;
    uint32_t uvOffset;

    // TODO: handle potentially 1+ uv sets
    uint32_t udimBase;
    uint32_t udimTableOffset;
};

struct DeviceScene
{
    Object *objects;
    Vec2f *float2s;
};

} // namespace ybi
