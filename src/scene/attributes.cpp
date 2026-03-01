#include "scene/attributes.h"

#include <cstring>

YBI_NAMESPACE_BEGIN

static Array<uint8_t> CopyBytes(MemoryView<uint8_t> view)
{
    Array<uint8_t> out(view.size());
    if (view.size() > 0)
    {
        memcpy(out.data(), view.data(), view.size());
    }
    return out;
}

static Array<int> CopyInts(MemoryView<int> view)
{
    Array<int> out(view.size());
    if (view.size() > 0)
    {
        memcpy(out.data(), view.data(), sizeof(int) * view.size());
    }
    return out;
}

Attribute::Attribute(Array<uint8_t> data,
                     AttributeType type,
                     PrimvarInterpolation interpolation,
                     const std::string &name)
    : name(name), data(std::move(data)), type(type), interpolation(interpolation)
{
}

Attribute::Attribute(Array<uint8_t> data,
                     Array<int> indices,
                     AttributeType type,
                     PrimvarInterpolation interpolation,
                     const std::string &name)
    : name(name), data(std::move(data)), indices(std::move(indices)), type(type),
      interpolation(interpolation)
{
}

Attribute::Attribute(MemoryView<uint8_t> data,
                     AttributeType type,
                     PrimvarInterpolation interpolation,
                     const std::string &name)
    : Attribute(CopyBytes(data), type, interpolation, name)
{
}

Attribute::Attribute(MemoryView<uint8_t> data,
                     MemoryView<int> indices,
                     AttributeType type,
                     PrimvarInterpolation interpolation,
                     const std::string &name)
    : Attribute(CopyBytes(data), CopyInts(indices), type, interpolation, name)
{
}

YBI_NAMESPACE_END
