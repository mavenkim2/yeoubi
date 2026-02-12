#include "scene/attributes.h"

YBI_NAMESPACE_BEGIN

Attribute::Attribute(MemoryView<uint8_t> data,
                     AttributeType type,
                     PrimvarInterpolation interpolation)
    : data(data), type(type), interpolation(interpolation)
{
}

YBI_NAMESPACE_END
