#include "core/assert.h"
#include "core/math.h"
#include "core/vec2f.h"
#include "core/vec3f.h"
#include "device/cuda_compat.h"
#include "render/bsdf/ltc_coefficients.h"
#include "render/sampling.h"
#include <cmath>

namespace ybi
{

YBI_GPU float ClampRandomNumber(float u)
{
    YBI_ASSERT(u >= 0.f);
    u = Min(u, FLOAT_BEFORE_ONE);
    return u;
}

} // namespace ybi
