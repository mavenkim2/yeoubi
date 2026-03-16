#pragma once

#include "util/vec3.h"

#include <cstdint>

namespace ybi
{

enum class UpAxis : uint8_t
{
    Y,
    Z,
};

YBI_DEVICE constexpr Vec3 UpAxisVector(UpAxis upAxis)
{
    return upAxis == UpAxis::Y ? Vec3(0.0f, 1.0f, 0.0f) : Vec3(0.0f, 0.0f, 1.0f);
}

} // namespace ybi
