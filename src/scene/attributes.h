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

struct Attribute
{
    MemoryView<uint8_t> data;
    AttributeType type;
    PrimvarInterpolation interpolation;

    Attribute(MemoryView<uint8_t> data, AttributeType type, PrimvarInterpolation interpolation);
};

YBI_NAMESPACE_END
