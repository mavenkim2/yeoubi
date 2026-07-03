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
};

struct DeviceScene
{
    Object *objects;
    Vec2f *float2s;
};

}
