#pragma once

#include "util/vec2.h"
#include "util/vec3.h"
#include "util/vec4.h"

#include <assert.h>
#include <math.h>
#include <stdint.h>

using vec2 = ybi::Vec2;
using vec3 = ybi::Vec3;
using vec4 = ybi::Vec4;
using uint = unsigned int;

template <typename T>
YBI_INTEGRATOR_HD constexpr T openpbr_ybi_min(const T a, const T b)
{
    return a < b ? a : b;
}

template <typename T>
YBI_INTEGRATOR_HD constexpr T openpbr_ybi_max(const T a, const T b)
{
    return a > b ? a : b;
}

YBI_INTEGRATOR_HD float openpbr_ybi_abs(const float x)
{
    return fabsf(x);
}

YBI_INTEGRATOR_HD vec2 openpbr_ybi_abs(const vec2 v)
{
    return vec2(fabsf(v.x), fabsf(v.y));
}

YBI_INTEGRATOR_HD vec3 openpbr_ybi_abs(const vec3 v)
{
    return vec3(fabsf(v.x), fabsf(v.y), fabsf(v.z));
}

YBI_INTEGRATOR_HD vec4 openpbr_ybi_abs(const vec4 v)
{
    return vec4(fabsf(v.x), fabsf(v.y), fabsf(v.z), fabsf(v.w));
}

YBI_INTEGRATOR_HD float openpbr_ybi_acos(const float x)
{
    return acosf(x);
}

YBI_INTEGRATOR_HD float openpbr_ybi_atan(const float x)
{
    return atanf(x);
}

YBI_INTEGRATOR_HD float openpbr_ybi_atan(const float y, const float x)
{
    return atan2f(y, x);
}

YBI_INTEGRATOR_HD float openpbr_ybi_cos(const float x)
{
    return cosf(x);
}

YBI_INTEGRATOR_HD float openpbr_ybi_exp(const float x)
{
    return expf(x);
}

YBI_INTEGRATOR_HD vec3 openpbr_ybi_exp(const vec3 v)
{
    return vec3(expf(v.x), expf(v.y), expf(v.z));
}

YBI_INTEGRATOR_HD float openpbr_ybi_floor(const float x)
{
    return floorf(x);
}

YBI_INTEGRATOR_HD float openpbr_ybi_log(const float x)
{
    return logf(x);
}

YBI_INTEGRATOR_HD vec3 openpbr_ybi_log(const vec3 v)
{
    return vec3(logf(v.x), logf(v.y), logf(v.z));
}

YBI_INTEGRATOR_HD float openpbr_ybi_pow(const float x, const float y)
{
    return powf(x, y);
}

YBI_INTEGRATOR_HD vec3 openpbr_ybi_pow(const vec3 a, const vec3 b)
{
    return vec3(powf(a.x, b.x), powf(a.y, b.y), powf(a.z, b.z));
}

YBI_INTEGRATOR_HD vec3 openpbr_ybi_pow(const vec3 a, const float b)
{
    return vec3(powf(a.x, b), powf(a.y, b), powf(a.z, b));
}

YBI_INTEGRATOR_HD float openpbr_ybi_sin(const float x)
{
    return sinf(x);
}

YBI_INTEGRATOR_HD float openpbr_ybi_sqrt(const float x)
{
    return sqrtf(x);
}

YBI_INTEGRATOR_HD vec3 openpbr_ybi_sqrt(const vec3 v)
{
    return vec3(sqrtf(v.x), sqrtf(v.y), sqrtf(v.z));
}

YBI_INTEGRATOR_HD vec2 openpbr_ybi_min(const vec2 a, const vec2 b)
{
    return vec2(openpbr_ybi_min(a.x, b.x), openpbr_ybi_min(a.y, b.y));
}

YBI_INTEGRATOR_HD vec3 openpbr_ybi_min(const vec3 a, const vec3 b)
{
    return vec3(openpbr_ybi_min(a.x, b.x), openpbr_ybi_min(a.y, b.y), openpbr_ybi_min(a.z, b.z));
}

YBI_INTEGRATOR_HD vec4 openpbr_ybi_min(const vec4 a, const vec4 b)
{
    return vec4(openpbr_ybi_min(a.x, b.x),
                openpbr_ybi_min(a.y, b.y),
                openpbr_ybi_min(a.z, b.z),
                openpbr_ybi_min(a.w, b.w));
}

YBI_INTEGRATOR_HD vec2 openpbr_ybi_max(const vec2 a, const vec2 b)
{
    return vec2(openpbr_ybi_max(a.x, b.x), openpbr_ybi_max(a.y, b.y));
}

YBI_INTEGRATOR_HD vec3 openpbr_ybi_max(const vec3 a, const vec3 b)
{
    return vec3(openpbr_ybi_max(a.x, b.x), openpbr_ybi_max(a.y, b.y), openpbr_ybi_max(a.z, b.z));
}

YBI_INTEGRATOR_HD vec4 openpbr_ybi_max(const vec4 a, const vec4 b)
{
    return vec4(openpbr_ybi_max(a.x, b.x),
                openpbr_ybi_max(a.y, b.y),
                openpbr_ybi_max(a.z, b.z),
                openpbr_ybi_max(a.w, b.w));
}

YBI_INTEGRATOR_HD float openpbr_ybi_clamp(const float x, const float lo, const float hi)
{
    return openpbr_ybi_max(lo, openpbr_ybi_min(x, hi));
}

YBI_INTEGRATOR_HD vec2 openpbr_ybi_clamp(const vec2 v, const vec2 lo, const vec2 hi)
{
    return openpbr_ybi_max(lo, openpbr_ybi_min(v, hi));
}

YBI_INTEGRATOR_HD vec3 openpbr_ybi_clamp(const vec3 v, const vec3 lo, const vec3 hi)
{
    return openpbr_ybi_max(lo, openpbr_ybi_min(v, hi));
}

YBI_INTEGRATOR_HD vec4 openpbr_ybi_clamp(const vec4 v, const vec4 lo, const vec4 hi)
{
    return openpbr_ybi_max(lo, openpbr_ybi_min(v, hi));
}

YBI_INTEGRATOR_HD vec2 openpbr_ybi_clamp(const vec2 v, const float lo, const float hi)
{
    return openpbr_ybi_clamp(v, vec2(lo), vec2(hi));
}

YBI_INTEGRATOR_HD vec3 openpbr_ybi_clamp(const vec3 v, const float lo, const float hi)
{
    return openpbr_ybi_clamp(v, vec3(lo), vec3(hi));
}

YBI_INTEGRATOR_HD vec4 openpbr_ybi_clamp(const vec4 v, const float lo, const float hi)
{
    return openpbr_ybi_clamp(v, vec4(lo), vec4(hi));
}

YBI_INTEGRATOR_HD float openpbr_ybi_mix(const float a, const float b, const float t)
{
    return a * (1.0f - t) + b * t;
}

YBI_INTEGRATOR_HD vec2 openpbr_ybi_mix(const vec2 a, const vec2 b, const float t)
{
    return a * (1.0f - t) + b * t;
}

YBI_INTEGRATOR_HD vec3 openpbr_ybi_mix(const vec3 a, const vec3 b, const float t)
{
    return a * (1.0f - t) + b * t;
}

YBI_INTEGRATOR_HD vec4 openpbr_ybi_mix(const vec4 a, const vec4 b, const float t)
{
    return a * (1.0f - t) + b * t;
}

YBI_INTEGRATOR_HD vec2 openpbr_ybi_mix(const vec2 a, const vec2 b, const vec2 t)
{
    return vec2(openpbr_ybi_mix(a.x, b.x, t.x), openpbr_ybi_mix(a.y, b.y, t.y));
}

YBI_INTEGRATOR_HD vec3 openpbr_ybi_mix(const vec3 a, const vec3 b, const vec3 t)
{
    return vec3(openpbr_ybi_mix(a.x, b.x, t.x),
                openpbr_ybi_mix(a.y, b.y, t.y),
                openpbr_ybi_mix(a.z, b.z, t.z));
}

YBI_INTEGRATOR_HD float openpbr_ybi_smoothstep(const float edge0, const float edge1, const float x)
{
    const float denom = openpbr_ybi_max(edge1 - edge0, 1.0e-8f);
    const float t = openpbr_ybi_clamp((x - edge0) / denom, 0.0f, 1.0f);
    return t * t * (3.0f - 2.0f * t);
}

YBI_INTEGRATOR_HD float openpbr_ybi_dot(const vec2 a, const vec2 b)
{
    return ybi::Dot(a, b);
}

YBI_INTEGRATOR_HD float openpbr_ybi_dot(const vec3 a, const vec3 b)
{
    return ybi::Dot(a, b);
}

YBI_INTEGRATOR_HD float openpbr_ybi_dot(const vec4 a, const vec4 b)
{
    return ybi::Dot(a, b);
}

YBI_INTEGRATOR_HD vec3 openpbr_ybi_cross(const vec3 a, const vec3 b)
{
    return ybi::Cross(a, b);
}

YBI_INTEGRATOR_HD float openpbr_ybi_length(const vec2 v)
{
    return ybi::Length(v);
}

YBI_INTEGRATOR_HD float openpbr_ybi_length(const vec3 v)
{
    return ybi::Length(v);
}

YBI_INTEGRATOR_HD float openpbr_ybi_length(const vec4 v)
{
    return ybi::Length(v);
}

YBI_INTEGRATOR_HD vec3 openpbr_ybi_normalize(const vec3 v)
{
    return ybi::Normalize(v);
}

YBI_INTEGRATOR_HD vec3 openpbr_ybi_reflect(const vec3 incident, const vec3 normal)
{
    return ybi::Reflect(incident, normal);
}

YBI_INTEGRATOR_HD vec2 openpbr_ybi_equal(const vec2 a, const vec2 b)
{
    return vec2(a.x == b.x ? 1.0f : 0.0f, a.y == b.y ? 1.0f : 0.0f);
}

YBI_INTEGRATOR_HD vec3 openpbr_ybi_equal(const vec3 a, const vec3 b)
{
    return vec3(a.x == b.x ? 1.0f : 0.0f, a.y == b.y ? 1.0f : 0.0f, a.z == b.z ? 1.0f : 0.0f);
}

YBI_INTEGRATOR_HD vec2 openpbr_ybi_greaterThan(const vec2 a, const vec2 b)
{
    return vec2(a.x > b.x ? 1.0f : 0.0f, a.y > b.y ? 1.0f : 0.0f);
}

YBI_INTEGRATOR_HD vec3 openpbr_ybi_greaterThan(const vec3 a, const vec3 b)
{
    return vec3(a.x > b.x ? 1.0f : 0.0f, a.y > b.y ? 1.0f : 0.0f, a.z > b.z ? 1.0f : 0.0f);
}

YBI_INTEGRATOR_HD vec3 openpbr_ybi_greaterThanEqual(const vec3 a, const vec3 b)
{
    return vec3(a.x >= b.x ? 1.0f : 0.0f, a.y >= b.y ? 1.0f : 0.0f, a.z >= b.z ? 1.0f : 0.0f);
}

YBI_INTEGRATOR_HD vec2 openpbr_ybi_notEqual(const vec2 a, const vec2 b)
{
    return vec2(a.x != b.x ? 1.0f : 0.0f, a.y != b.y ? 1.0f : 0.0f);
}

YBI_INTEGRATOR_HD vec3 openpbr_ybi_notEqual(const vec3 a, const vec3 b)
{
    return vec3(a.x != b.x ? 1.0f : 0.0f, a.y != b.y ? 1.0f : 0.0f, a.z != b.z ? 1.0f : 0.0f);
}

YBI_INTEGRATOR_HD vec4 openpbr_ybi_notEqual(const vec4 a, const vec4 b)
{
    return vec4(a.x != b.x ? 1.0f : 0.0f,
                a.y != b.y ? 1.0f : 0.0f,
                a.z != b.z ? 1.0f : 0.0f,
                a.w != b.w ? 1.0f : 0.0f);
}

YBI_INTEGRATOR_HD bool openpbr_ybi_all(const vec2 v)
{
    return v.x != 0.0f && v.y != 0.0f;
}

YBI_INTEGRATOR_HD bool openpbr_ybi_all(const vec3 v)
{
    return v.x != 0.0f && v.y != 0.0f && v.z != 0.0f;
}

YBI_INTEGRATOR_HD bool openpbr_ybi_any(const vec2 v)
{
    return v.x != 0.0f || v.y != 0.0f;
}

YBI_INTEGRATOR_HD bool openpbr_ybi_any(const vec3 v)
{
    return v.x != 0.0f || v.y != 0.0f || v.z != 0.0f;
}

YBI_INTEGRATOR_HD bool openpbr_ybi_any(const vec4 v)
{
    return v.x != 0.0f || v.y != 0.0f || v.z != 0.0f || v.w != 0.0f;
}

YBI_INTEGRATOR_HD float openpbr_ybi_saturate(const float x)
{
    return openpbr_ybi_clamp(x, 0.0f, 1.0f);
}

YBI_INTEGRATOR_HD vec3 openpbr_ybi_saturate(const vec3 v)
{
    return openpbr_ybi_clamp(v, 0.0f, 1.0f);
}

YBI_INTEGRATOR_HD vec2 openpbr_swizzle_xy(const vec3 v)
{
    return vec2(v.x, v.y);
}

YBI_INTEGRATOR_HD vec2 openpbr_swizzle_xy(const vec4 v)
{
    return vec2(v.x, v.y);
}

YBI_INTEGRATOR_HD vec3 openpbr_swizzle_xyz(const vec4 v)
{
    return vec3(v.x, v.y, v.z);
}

#define OPENPBR_USE_CUSTOM_INTEROP 1

#define ADDRESS_SPACE_THREAD
#define OUT(type) type&
#define INOUT(type) type&
#define CONST_REF(type) const type&
#define CONSTEXPR_LOCAL constexpr
#define CONSTEXPR_GLOBAL static inline constexpr
#define GENERAL_CONSTEXPR_FUNCTION static constexpr
#define LIMITED_CONSTEXPR_FUNCTION static constexpr
#if defined(__CUDACC__)
#define INLINE_FUNCTION __host__ __device__ inline
#else
#define INLINE_FUNCTION inline
#endif
#define SWIZZLE(v, suffix) openpbr_swizzle_##suffix(v)

#define MAKE_STRUCT_1(type, arg1) type{arg1}
#define MAKE_STRUCT_2(type, arg1, arg2) type{arg1, arg2}
#define MAKE_STRUCT_3(type, arg1, arg2, arg3) type{arg1, arg2, arg3}
#define MAKE_STRUCT_4(type, arg1, arg2, arg3, arg4) type{arg1, arg2, arg3, arg4}
#define MAKE_STRUCT_5(type, arg1, arg2, arg3, arg4, arg5) type{arg1, arg2, arg3, arg4, arg5}
#define MAKE_STRUCT_6(type, arg1, arg2, arg3, arg4, arg5, arg6) type{arg1, arg2, arg3, arg4, arg5, arg6}
#define MAKE_STRUCT_7(type, arg1, arg2, arg3, arg4, arg5, arg6, arg7) type{arg1, arg2, arg3, arg4, arg5, arg6, arg7}
#define MAKE_STRUCT_8(type, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8) type{arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8}
#define MAKE_STRUCT_9(type, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9) type{arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9}
#define MAKE_STRUCT_10(type, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10) type{arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10}
#define MAKE_STRUCT_11(type, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10, arg11) type{arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10, arg11}
#define MAKE_STRUCT_12(type, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10, arg11, arg12) type{arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10, arg11, arg12}
#define MAKE_STRUCT_13(type, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10, arg11, arg12, arg13) type{arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10, arg11, arg12, arg13}
#define MAKE_STRUCT_14(type, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10, arg11, arg12, arg13, arg14) type{arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10, arg11, arg12, arg13, arg14}
#define MAKE_STRUCT_15(type, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10, arg11, arg12, arg13, arg14, arg15) type{arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10, arg11, arg12, arg13, arg14, arg15}

#define DECLARE_SPECIALIZATION_CONSTANT(constant_id_number, name, default_value) CONSTEXPR_GLOBAL bool name = default_value
#define GET_SPECIALIZATION_CONSTANT(name) name

#if defined(__CUDACC__)
#define ASSERT(expr, message) ((void)0)
#define ASSERT_UNREACHABLE(message) ((void)0)
#else
#define ASSERT(expr, message) assert((expr) && (message))
#define ASSERT_UNREACHABLE(message) assert(false && (message))
#endif
#define STATIC_ASSERT(expr, message) static_assert(expr, message)

#define abs openpbr_ybi_abs
#define acos openpbr_ybi_acos
#define atan openpbr_ybi_atan
#define clamp openpbr_ybi_clamp
#define cos openpbr_ybi_cos
#define exp openpbr_ybi_exp
#define floor openpbr_ybi_floor
#define log openpbr_ybi_log
#define max openpbr_ybi_max
#define min openpbr_ybi_min
#define mix openpbr_ybi_mix
#define pow openpbr_ybi_pow
#define sin openpbr_ybi_sin
#define smoothstep openpbr_ybi_smoothstep
#define sqrt openpbr_ybi_sqrt
#define cross openpbr_ybi_cross
#define dot openpbr_ybi_dot
#define length openpbr_ybi_length
#define normalize openpbr_ybi_normalize
#define reflect openpbr_ybi_reflect
#define equal openpbr_ybi_equal
#define greaterThan openpbr_ybi_greaterThan
#define greaterThanEqual openpbr_ybi_greaterThanEqual
#define notEqual openpbr_ybi_notEqual
#define all openpbr_ybi_all
#define any openpbr_ybi_any
#define saturate openpbr_ybi_saturate
