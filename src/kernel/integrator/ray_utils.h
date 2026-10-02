#pragma once

#include "core/bit.h"
#include "core/vec3f.h"
#include "core/vec3i.h"

namespace ybi
{

YBI_CPU_GPU Vec3f OffsetRayOrigin(const Vec3f &p, const Vec3f &gn, bool printDebug = false)
{
    static const float intScale = 256.f;
    static const float floatScale = 1.f / 65536;
    static const float origin = 1.f / 32;

    Vec3f of_f = intScale * gn;
    Vec3i of_i(int(of_f.x), int(of_f.y), int(of_f.z));

    Vec3f p_i = Vec3f(IntAsFloat(FloatAsInt(p.x) + ((p.x < 0) ? -of_i.x : of_i.x)),
                      IntAsFloat(FloatAsInt(p.y) + ((p.y < 0) ? -of_i.y : of_i.y)),
                      IntAsFloat(FloatAsInt(p.z) + ((p.z < 0) ? -of_i.z : of_i.z)));

    return Vec3f(Abs(p.x) < origin ? p.x + floatScale * gn.x : p_i.x,
                 Abs(p.y) < origin ? p.y + floatScale * gn.y : p_i.y,
                 Abs(p.z) < origin ? p.z + floatScale * gn.z : p_i.z);
}

} // namespace ybi
