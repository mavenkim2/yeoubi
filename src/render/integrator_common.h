#pragma once

#include "util/math_constants.h"
#include "util/vec2.h"
#include "util/vec3.h"
#include "util/vec4.h"
#include <cmath>
#include <cstdint>

namespace ybi
{
namespace render
{
namespace integrator
{

using ybi::BuildOrthonormalBasis;
using ybi::Clamp;
using ybi::Cross;
using ybi::Dot;
using ybi::FaceForward;
using ybi::Vec2;
using ybi::Vec3;
using ybi::Vec4;
using ybi::Length;
using ybi::LengthSquared;
using ybi::Lerp;
using ybi::Luminance;
using ybi::MaxComponent;
using ybi::Normalize;
using ybi::Reflect;

static constexpr int kSemanticDiffuse = 0;
static constexpr int kSemanticRoughness = 1;
static constexpr int kSemanticMetallic = 2;
static constexpr int kSemanticOcclusion = 3;
static constexpr int kSemanticNormal = 4;
static constexpr int kSemanticIor = 5;
static constexpr int kSemanticEmissive = 6;
static constexpr int kSemanticOpacity = 7;
static constexpr int kSemanticSpecularColor = 8;
static constexpr int kSemanticClearcoat = 9;
static constexpr int kSemanticClearcoatRoughness = 10;

static constexpr int kWrapModeUnknown = 0;
static constexpr int kWrapModeRepeat = 1;
static constexpr int kWrapModeClamp = 2;
static constexpr int kWrapModeMirror = 3;
static constexpr int kWrapModeBlack = 4;
static constexpr int kWrapModeUseMetadata = 5;

static constexpr unsigned int kVirtualTexturePageTypeStream = 1u;
static constexpr unsigned int kVirtualTexturePageTypeTail = 2u;

struct UInt2
{
    unsigned int x;
    unsigned int y;
};

YBI_DEVICE bool ShouldRenderLaunchPixel(
    int singlePixelEnabled, int singlePixelX, int singlePixelY, unsigned int x, unsigned int y)
{
    if (singlePixelEnabled == 0)
    {
        return true;
    }
    if (singlePixelX < 0 || singlePixelY < 0)
    {
        return false;
    }
    return x == static_cast<unsigned int>(singlePixelX) &&
           y == static_cast<unsigned int>(singlePixelY);
}

struct UV2
{
    float x;
    float y;
};

struct HitInfo
{
    int instanceId = -1;
    int primitiveIndex = -1;
    Vec3 barycentrics = Vec3(0.0f, 0.0f, 0.0f);
    Vec3 rayOrigin = Vec3(0.0f, 0.0f, 0.0f);
    Vec3 rayDir = Vec3(0.0f, 0.0f, 1.0f);
    float t = 0.0f;
    Vec3 geomNormal = Vec3(0.0f, 0.0f, 1.0f);
    Vec3 shadingNormal = Vec3(0.0f, 0.0f, 1.0f);
    Vec3 worldTri0 = Vec3(0.0f, 0.0f, 0.0f);
    Vec3 worldTri1 = Vec3(0.0f, 0.0f, 0.0f);
    Vec3 worldTri2 = Vec3(0.0f, 0.0f, 0.0f);
    Vec3 dPds = Vec3(0.0f, 0.0f, 0.0f);
    Vec3 dPdt = Vec3(0.0f, 0.0f, 0.0f);
    Vec3 dpdx = Vec3(0.0f, 0.0f, 0.0f);
    Vec3 dpdy = Vec3(0.0f, 0.0f, 0.0f);
    float dSdx = 0.0f;
    float dTdx = 0.0f;
    float dSdy = 0.0f;
    float dTdy = 0.0f;
    bool hasBarycentrics = false;
    bool hasGeomNormal = false;
    bool hasShadingNormal = false;
    bool hasWorldTriangle = false;
    bool hasSurfacePartials = false;
    bool hasPositionDifferentials = false;
    bool hasTextureDifferentials = false;
};

YBI_DEVICE float Clamp01(float v)
{
    if (v < 0.0f)
    {
        return 0.0f;
    }
    if (v > 1.0f)
    {
        return 1.0f;
    }
    return v;
}

YBI_DEVICE int ClampInt(int v, int lo, int hi)
{
    return v < lo ? lo : (v > hi ? hi : v);
}

YBI_DEVICE int MaxInt(int a, int b)
{
    return a > b ? a : b;
}

YBI_DEVICE int32_t FloatAsInt(float value)
{
    union
    {
        float f;
        int32_t i;
    } bits = {value};
    return bits.i;
}

YBI_DEVICE float IntAsFloat(int32_t value)
{
    union
    {
        int32_t i;
        float f;
    } bits = {value};
    return bits.f;
}

YBI_DEVICE Vec3 OffsetRayOrigin(const Vec3 &point, const Vec3 &normal)
{
    constexpr float kOrigin = 1.0f / 32.0f;
    constexpr float kFloatScale = 1.0f / 65536.0f;
    constexpr float kIntScale = 256.0f;

    const int32_t ofx = static_cast<int32_t>(kIntScale * normal.x);
    const int32_t ofy = static_cast<int32_t>(kIntScale * normal.y);
    const int32_t ofz = static_cast<int32_t>(kIntScale * normal.z);

    const Vec3 offsetPoint(
        IntAsFloat(FloatAsInt(point.x) + (point.x < 0.0f ? -ofx : ofx)),
        IntAsFloat(FloatAsInt(point.y) + (point.y < 0.0f ? -ofy : ofy)),
        IntAsFloat(FloatAsInt(point.z) + (point.z < 0.0f ? -ofz : ofz)));

    return Vec3(fabsf(point.x) < kOrigin ? point.x + kFloatScale * normal.x : offsetPoint.x,
                fabsf(point.y) < kOrigin ? point.y + kFloatScale * normal.y : offsetPoint.y,
                fabsf(point.z) < kOrigin ? point.z + kFloatScale * normal.z : offsetPoint.z);
}

YBI_DEVICE Vec3 OffsetRayOrigin(const Vec3 &point,
                                       const Vec3 &geomNormal,
                                       const Vec3 &rayDirection)
{
    const Vec3 offsetNormal =
        Dot(geomNormal, rayDirection) >= 0.0f ? geomNormal : -geomNormal;
    return OffsetRayOrigin(point, offsetNormal);
}

YBI_DEVICE unsigned int Hash32(unsigned int x)
{
    x ^= x >> 16;
    x *= 0x7feb352du;
    x ^= x >> 15;
    x *= 0x846ca68bu;
    x ^= x >> 16;
    return x;
}

YBI_DEVICE float Random01(unsigned int &state)
{
    state = Hash32(state + 0x9e3779b9u);
    return float(state & 0x00ffffffu) / float(0x01000000u);
}

YBI_DEVICE Vec3 SampleCosineHemisphere(float u1, float u2)
{
    const float r = sqrtf(u1 < 0.0f ? 0.0f : u1);
    const float phi = ybi::kTwoPi * u2;
    const float x = r * cosf(phi);
    const float y = r * sinf(phi);
    const float z = sqrtf((1.0f - u1) < 0.0f ? 0.0f : (1.0f - u1));
    return Vec3(x, y, z);
}

YBI_DEVICE Vec3 SkyColor(const Vec3 &direction)
{
    const float t = 0.5f * (direction.y + 1.0f);
    const Vec3 top(0.7f, 0.8f, 1.0f);
    const Vec3 bottom(0.2f, 0.25f, 0.35f);
    return Vec3((1.0f - t) * top.x + t * bottom.x,
                (1.0f - t) * top.y + t * bottom.y,
                (1.0f - t) * top.z + t * bottom.z);
}

YBI_DEVICE float ApplyWrapMode(float uv, int wrapMode, bool &outBlack)
{
    outBlack = false;
    switch (wrapMode)
    {
        case kWrapModeClamp:
            return Clamp01(uv);
        case kWrapModeMirror:
        {
            float t = fmodf(uv, 2.0f);
            if (t < 0.0f)
            {
                t += 2.0f;
            }
            const float mirrored = (t <= 1.0f) ? t : (2.0f - t);
            return Clamp01(mirrored);
        }
        case kWrapModeBlack:
            if (uv < 0.0f || uv > 1.0f)
            {
                outBlack = true;
                return 0.0f;
            }
            return uv;
        case kWrapModeRepeat:
        case kWrapModeUseMetadata:
        case kWrapModeUnknown:
        default:
            return uv - floorf(uv);
    }
}

} // namespace integrator
} // namespace render
} // namespace ybi
