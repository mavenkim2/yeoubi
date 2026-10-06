#pragma once

#include "core/vec3f.h"

namespace ybi
{

// Each coordinate holds the same coordinate from several independent samples.
template <typename Float>
struct LaneVec3
{
    Float x, y, z;

    LaneVec3() = default;
    YBI_GPU explicit LaneVec3(Float value) : x(value), y(value), z(value) {}
    YBI_GPU LaneVec3(Float x, Float y, Float z) : x(x), y(y), z(z) {}
    YBI_GPU explicit LaneVec3(const Vec3f &value) : x(value.x), y(value.y), z(value.z) {}

    YBI_GPU LaneVec3 &operator+=(const LaneVec3 &other)
    {
        x += other.x;
        y += other.y;
        z += other.z;
        return *this;
    }
};

template <typename Float>
YBI_GPU LaneVec3<Float> operator*(const LaneVec3<Float> &value, Float scale)
{
    return LaneVec3<Float>(value.x * scale, value.y * scale, value.z * scale);
}

template <typename Float>
struct Vector3Traits;

template <>
struct Vector3Traits<float>
{
    using Type = Vec3f;
};

template <typename Float>
using Vector3 = typename Vector3Traits<Float>::Type;

YBI_CPU_GPU float Select(bool mask, float ifTrue, float ifFalse)
{
    return mask ? ifTrue : ifFalse;
}

YBI_CPU_GPU bool All(bool mask)
{
    return mask;
}

} // namespace ybi
