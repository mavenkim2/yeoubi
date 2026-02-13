#pragma once

#include "util/memory_view.h"
#include <cstdint>

YBI_NAMESPACE_BEGIN

enum class PrimvarInterpolation
{
    Constant,
    Uniform,
    Vertex,
    Varying,
    FaceVarying,
    Unknown
};

enum class AttributeType
{
    Float,
    Float2,
    Float3,
    Float4,
    Int,
    Int8,
    Bool,
    String,
    Quaternion,
    Unknown
};

inline bool AttributeTypeIsFloat(AttributeType type)
{
    return type == AttributeType::Float || type == AttributeType::Float2 ||
           type == AttributeType::Float3 || type == AttributeType::Float4;
}

inline size_t AttributeTypeGetSize(AttributeType type)
{
    switch (type)
    {
        case AttributeType::Float:
            return 4;
        case AttributeType::Float2:
            return 8;
        case AttributeType::Float3:
            return 12;
        case AttributeType::Float4:
            return 16;
        case AttributeType::Int:
            return 4;
        case AttributeType::Int8:
            return 1;
        case AttributeType::Bool:
            return 1;
        case AttributeType::String:
            return 0;
        case AttributeType::Quaternion:
            return 16;
        default:
            return 0;
    }
}

struct Attribute
{
    MemoryView<uint8_t> data;
    MemoryView<int> indices;
    AttributeType type;
    PrimvarInterpolation interpolation;

    Attribute(MemoryView<uint8_t> data, AttributeType type, PrimvarInterpolation interpolation);
    Attribute(MemoryView<uint8_t> data,
              MemoryView<int> indices,
              AttributeType type,
              PrimvarInterpolation interpolation);
};

YBI_NAMESPACE_END
