#pragma once

#include "core/assert.h"
#include "core/vec3f.h"
#include "device/cuda_compat.h"
#include <cmath>

namespace ybi
{

struct OrthonormalBasis
{
    // NOTE: column vectors
    Vec3f x;
    Vec3f y;
    Vec3f z;

    YBI_CPU_GPU OrthonormalBasis(const Vec3f &x, const Vec3f &y, const Vec3f &z) : x(x), y(y), z(z)
    {
    }
    YBI_CPU_GPU OrthonormalBasis(const Vec3f &x, const Vec3f &z)
        : x(Normalize(x)), y(Normalize(Cross(Normalize(z), Normalize(x)))), z(Normalize(z))
    {
    }

    YBI_CPU_GPU Vec3f LocalToWorld(const Vec3f &v) const
    {
        Vec3f result = x * v.x + y * v.y + z * v.z;
        return result;
    }

    YBI_CPU_GPU Vec3f WorldToLocal(const Vec3f &v) const
    {
        Vec3f result(Dot(v, x), Dot(v, y), Dot(v, z));
        return result;
    }
};

YBI_CPU_GPU void BuildOrthonormalBasis(const Vec3f &n, Vec3f &t0, Vec3f &t1)
{
    const float s = n.z >= 0.f ? 1.f : -1.f;
    const float a = -1.f / (s + n.z);
    const float b = n.x * n.y * a;

    t0 = Vec3f(fmaf(n.x * a, n.x * s, 1.0f), s * b, -s * n.x);
    t1 = Vec3f(b, fmaf(n.y * a, n.y, s), -n.y);
}

} // namespace ybi
