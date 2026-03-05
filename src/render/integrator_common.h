#pragma once

#include <cmath>
#include <cstdint>

#if defined(__CUDACC__)
#define YBI_INTEGRATOR_HD __host__ __device__ inline
#else
#define YBI_INTEGRATOR_HD inline
#endif

namespace ybi
{
namespace render
{
namespace integrator
{

static constexpr int kSemanticDiffuse = 0;
static constexpr int kSemanticRoughness = 1;
static constexpr int kSemanticMetallic = 2;
static constexpr int kSemanticOcclusion = 3;
static constexpr int kSemanticNormal = 4;
static constexpr int kSemanticIor = 5;
static constexpr int kSemanticEmissive = 6;
static constexpr int kSemanticOpacity = 7;

static constexpr int kWrapModeUnknown = 0;
static constexpr int kWrapModeRepeat = 1;
static constexpr int kWrapModeClamp = 2;
static constexpr int kWrapModeMirror = 3;
static constexpr int kWrapModeBlack = 4;
static constexpr int kWrapModeUseMetadata = 5;

static constexpr unsigned long long kVirtualTextureEmptyKey = ~0ull;
static constexpr unsigned int kVirtualTexturePageTypeInvalid = 0u;
static constexpr unsigned int kVirtualTexturePageTypeStream = 1u;
static constexpr unsigned int kVirtualTexturePageTypeTail = 2u;

struct UInt2
{
    unsigned int x;
    unsigned int y;
};

struct Vec2
{
    float x;
    float y;
};

struct Vec3
{
    float x;
    float y;
    float z;
};

struct Vec4
{
    float x;
    float y;
    float z;
    float w;
};

struct UV2
{
    float x;
    float y;
};

struct HitInfo
{
    int instanceId = -1;
    int primitiveIndex = -1;
    Vec3 barycentrics = {0.0f, 0.0f, 0.0f};
    Vec3 rayOrigin = {0.0f, 0.0f, 0.0f};
    Vec3 rayDir = {0.0f, 0.0f, 1.0f};
    float t = 0.0f;
    Vec3 geomNormal = {0.0f, 0.0f, 1.0f};
    bool hasBarycentrics = false;
    bool hasGeomNormal = false;
};

YBI_INTEGRATOR_HD Vec3 MakeVec3(float x, float y, float z)
{
    return {x, y, z};
}

YBI_INTEGRATOR_HD Vec4 MakeVec4(float x, float y, float z, float w)
{
    return {x, y, z, w};
}

YBI_INTEGRATOR_HD Vec3 Add(const Vec3 &a, const Vec3 &b)
{
    return {a.x + b.x, a.y + b.y, a.z + b.z};
}

YBI_INTEGRATOR_HD Vec3 Sub(const Vec3 &a, const Vec3 &b)
{
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}

YBI_INTEGRATOR_HD Vec3 Mul(const Vec3 &a, float b)
{
    return {a.x * b, a.y * b, a.z * b};
}

YBI_INTEGRATOR_HD float Dot(const Vec3 &a, const Vec3 &b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

YBI_INTEGRATOR_HD Vec3 Cross(const Vec3 &a, const Vec3 &b)
{
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}

YBI_INTEGRATOR_HD float Length(const Vec3 &a)
{
    return sqrtf(Dot(a, a));
}

YBI_INTEGRATOR_HD Vec3 Normalize(const Vec3 &v)
{
    const float len = Length(v);
    if (len <= 1e-20f)
    {
        return {0.0f, 0.0f, 1.0f};
    }
    const float inv = 1.0f / len;
    return {v.x * inv, v.y * inv, v.z * inv};
}

YBI_INTEGRATOR_HD Vec3 FaceForward(const Vec3 &normal, const Vec3 &referenceDirection)
{
    const float d = Dot(normal, referenceDirection);
    if (d < 0.0f)
    {
        return normal;
    }
    return {-normal.x, -normal.y, -normal.z};
}

YBI_INTEGRATOR_HD float Clamp01(float v)
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

YBI_INTEGRATOR_HD int ClampInt(int v, int lo, int hi)
{
    return v < lo ? lo : (v > hi ? hi : v);
}

YBI_INTEGRATOR_HD int MaxInt(int a, int b)
{
    return a > b ? a : b;
}

YBI_INTEGRATOR_HD unsigned int Hash32(unsigned int x)
{
    x ^= x >> 16;
    x *= 0x7feb352du;
    x ^= x >> 15;
    x *= 0x846ca68bu;
    x ^= x >> 16;
    return x;
}

YBI_INTEGRATOR_HD float Random01(unsigned int &state)
{
    state = Hash32(state + 0x9e3779b9u);
    return float(state & 0x00ffffffu) / float(0x01000000u);
}

YBI_INTEGRATOR_HD void BuildOrthonormalBasis(const Vec3 &n, Vec3 &t, Vec3 &b)
{
    const Vec3 up = fabsf(n.z) < 0.999f ? MakeVec3(0.0f, 0.0f, 1.0f)
                                        : MakeVec3(0.0f, 1.0f, 0.0f);
    t = Normalize(Cross(up, n));
    b = Normalize(Cross(n, t));
}

YBI_INTEGRATOR_HD Vec3 SampleCosineHemisphere(float u1, float u2)
{
    const float r = sqrtf(u1 < 0.0f ? 0.0f : u1);
    const float phi = 6.28318530718f * u2;
    const float x = r * cosf(phi);
    const float y = r * sinf(phi);
    const float z = sqrtf((1.0f - u1) < 0.0f ? 0.0f : (1.0f - u1));
    return MakeVec3(x, y, z);
}

YBI_INTEGRATOR_HD Vec3 SkyColor(const Vec3 &direction)
{
    const float t = 0.5f * (direction.y + 1.0f);
    const Vec3 top = MakeVec3(0.7f, 0.8f, 1.0f);
    const Vec3 bottom = MakeVec3(0.2f, 0.25f, 0.35f);
    return MakeVec3((1.0f - t) * top.x + t * bottom.x,
                    (1.0f - t) * top.y + t * bottom.y,
                    (1.0f - t) * top.z + t * bottom.z);
}

YBI_INTEGRATOR_HD float ApplyWrapMode(float uv, int wrapMode, bool &outBlack)
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
