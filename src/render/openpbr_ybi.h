#pragma once

#include "render/openpbr_ybi_interop.h"

#include "openpbr.h"

#undef abs
#undef acos
#undef atan
#undef clamp
#undef cos
#undef exp
#undef floor
#undef log
#undef max
#undef min
#undef mix
#undef pow
#undef sin
#undef smoothstep
#undef sqrt
#undef cross
#undef dot
#undef length
#undef normalize
#undef reflect
#undef equal
#undef greaterThan
#undef greaterThanEqual
#undef notEqual
#undef all
#undef any
#undef saturate

namespace ybi
{
namespace render
{
namespace integrator
{

YBI_INTEGRATOR_HD Vec3 OpenPbrDefaultRgbWavelengthsNm()
{
    return Vec3(700.0f, 546.1f, 435.8f);
}

YBI_INTEGRATOR_HD Vec3 OpenPbrCombinedWeight(const OpenPBR_DiffuseSpecular &weight)
{
    return Vec3(weight.diffuse.x + weight.specular.x,
                weight.diffuse.y + weight.specular.y,
                weight.diffuse.z + weight.specular.z);
}

YBI_INTEGRATOR_HD bool OpenPbrSampleIsDelta(const OpenPBR_BsdfLobeType sampledType)
{
    return (sampledType & OpenPBR_BsdfLobeTypeSpecular) != 0u;
}

} // namespace integrator
} // namespace render
} // namespace ybi

#undef OPENPBR_USE_CUSTOM_INTEROP
#undef ADDRESS_SPACE_THREAD
#undef OUT
#undef INOUT
#undef CONST_REF
#undef CONSTEXPR_LOCAL
#undef CONSTEXPR_GLOBAL
#undef GENERAL_CONSTEXPR_FUNCTION
#undef LIMITED_CONSTEXPR_FUNCTION
#undef INLINE_FUNCTION
#undef SWIZZLE
#undef MAKE_STRUCT_1
#undef MAKE_STRUCT_2
#undef MAKE_STRUCT_3
#undef MAKE_STRUCT_4
#undef MAKE_STRUCT_5
#undef MAKE_STRUCT_6
#undef MAKE_STRUCT_7
#undef MAKE_STRUCT_8
#undef MAKE_STRUCT_9
#undef MAKE_STRUCT_10
#undef MAKE_STRUCT_11
#undef MAKE_STRUCT_12
#undef MAKE_STRUCT_13
#undef MAKE_STRUCT_14
#undef MAKE_STRUCT_15
#undef DECLARE_SPECIALIZATION_CONSTANT
#undef GET_SPECIALIZATION_CONSTANT
#undef ASSERT
#undef ASSERT_UNREACHABLE
#undef STATIC_ASSERT
