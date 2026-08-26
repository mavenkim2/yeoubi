#include "core/vec2f.h"
#include "core/vec3f.h"
#include "device/cuda_compat.h"
#include "render/sampling.h"

namespace ybi
{

YBI_CPU_GPU Vec3f SampleSheenLTC(Vec2f u)
{
    Vec3f w = SampleCosineHemisphere(u);

    float aInv, bInv; // <-- from ltc table, inverse
    float a = 1.f / aInv;

    Vec3f result(a * w.x - bInv * a * w.z, a * w.y, w.z);
    return result;
}

YBI_CPU_GPU Vec3f ComeOn() {}

} // namespace ybi
