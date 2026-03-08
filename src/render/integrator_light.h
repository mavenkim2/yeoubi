#pragma once

#include "render/integrator_common.h"
#include "render/launch_params.h"
#include "scene/material_light.h"

namespace ybi
{
namespace render
{
namespace integrator
{

struct DirectLightSample
{
    int lightIndex = -1;
    Vec3 wi = MakeVec3(0.0f, 0.0f, 1.0f);
    Vec3 radiance = MakeVec3(0.0f, 0.0f, 0.0f);
    Vec3 lightPoint = MakeVec3(0.0f, 0.0f, 0.0f);
    Vec3 lightNormal = MakeVec3(0.0f, 0.0f, 1.0f);
    float distance = 0.0f;
    float pdf = 0.0f;
    bool isDeltaLight = false;
};

struct LightRayHit
{
    int lightIndex = -1;
    Vec3 radiance = MakeVec3(0.0f, 0.0f, 0.0f);
    float distance = 0.0f;
    float pdf = 0.0f;
    bool isDeltaLight = false;
};

YBI_INTEGRATOR_HD float MaxF(float a, float b)
{
    return a > b ? a : b;
}

YBI_INTEGRATOR_HD float MinF(float a, float b)
{
    return a < b ? a : b;
}

YBI_INTEGRATOR_HD Vec3 ToVec3(const ybi::float3 &value)
{
    return MakeVec3(value.x, value.y, value.z);
}

YBI_INTEGRATOR_HD const PackedLight *GetPackedLights(const LaunchParams &params);
YBI_INTEGRATOR_HD const int *GetLightShadowExcludeRefs(const LaunchParams &params);
YBI_INTEGRATOR_HD bool IsRefExcludedFromLightShadow(const LaunchParams &params,
                                                    int lightIndex,
                                                    int refIndex);
YBI_INTEGRATOR_HD bool LightHasShadowExcludes(const LaunchParams &params, int lightIndex);

YBI_INTEGRATOR_HD Vec3 LightEmission(const PackedLight &light)
{
    return Mul(ToVec3(light.color), light.emissionScale);
}

YBI_INTEGRATOR_HD bool IsFiniteAreaLight(const PackedLight &light)
{
    return light.type == static_cast<unsigned int>(LightType::Rect) ||
           light.type == static_cast<unsigned int>(LightType::Disk) ||
           light.type == static_cast<unsigned int>(LightType::Sphere) ||
           light.type == static_cast<unsigned int>(LightType::Cylinder);
}

YBI_INTEGRATOR_HD void GetLightBasis(const PackedLight &light, Vec3 *outT, Vec3 *outB, Vec3 *outN)
{
    Vec3 n = Normalize(ToVec3(light.direction));
    Vec3 t = Normalize(ToVec3(light.tangent));
    Vec3 b = Normalize(ToVec3(light.bitangent));
    if (fabsf(Dot(n, t)) > 0.999f || fabsf(Dot(n, b)) > 0.999f)
    {
        BuildOrthonormalBasis(n, t, b);
    }
    else
    {
        b = Normalize(Cross(n, t));
        t = Normalize(Cross(b, n));
    }
    if (outT)
    {
        *outT = t;
    }
    if (outB)
    {
        *outB = b;
    }
    if (outN)
    {
        *outN = n;
    }
}

YBI_INTEGRATOR_HD Vec3 EnvironmentRadiance(const LaunchParams &params, const Vec3 &direction)
{
    Vec3 env = MakeVec3(0.0f, 0.0f, 0.0f);
    bool hasDome = false;
    const PackedLight *lights = GetPackedLights(params);
    if (lights)
    {
        for (int i = 0; i < params.lightCount; ++i)
        {
            const PackedLight &light = lights[i];
            if (light.type != static_cast<unsigned int>(LightType::Dome))
            {
                continue;
            }
            env = Add(env, LightEmission(light));
            hasDome = true;
        }
    }
    if (!hasDome)
    {
        return SkyColor(direction);
    }
    return env;
}

YBI_INTEGRATOR_HD int CountDirectLights(const LaunchParams &params)
{
    const PackedLight *lights = GetPackedLights(params);
    if (!lights)
    {
        return 0;
    }

    int count = 0;
    for (int i = 0; i < params.lightCount; ++i)
    {
        if (lights[i].type != static_cast<unsigned int>(LightType::Dome))
        {
            count++;
        }
    }
    return count;
}

YBI_INTEGRATOR_HD float DirectLightPickPdf(const LaunchParams &params)
{
    const int count = CountDirectLights(params);
    return count > 0 ? 1.0f / float(count) : 0.0f;
}

YBI_INTEGRATOR_HD const PackedLight *GetPackedLights(const LaunchParams &params)
{
    if (params.lights == 0ull || params.lightCount <= 0)
    {
        return nullptr;
    }
    return reinterpret_cast<const PackedLight *>(params.lights);
}

YBI_INTEGRATOR_HD const int *GetLightShadowExcludeRefs(const LaunchParams &params)
{
    if (params.lightShadowExcludeRefs == 0ull || params.lightShadowExcludeRefCount <= 0)
    {
        return nullptr;
    }
    return reinterpret_cast<const int *>(params.lightShadowExcludeRefs);
}

YBI_INTEGRATOR_HD bool IsRefExcludedFromLightShadow(const LaunchParams &params,
                                                    int lightIndex,
                                                    int refIndex)
{
    const PackedLight *lights = GetPackedLights(params);
    const int *excludeRefs = GetLightShadowExcludeRefs(params);
    if (!lights || !excludeRefs || lightIndex < 0 || lightIndex >= params.lightCount || refIndex < 0)
    {
        return false;
    }

    const PackedLight &light = lights[lightIndex];
    if (light.shadowExcludeCount == 0u)
    {
        return false;
    }

    const uint32_t begin = light.shadowExcludeOffset;
    const uint32_t end = begin + light.shadowExcludeCount;
    if (end > static_cast<uint32_t>(params.lightShadowExcludeRefCount))
    {
        return false;
    }

    for (uint32_t i = begin; i < end; ++i)
    {
        if (excludeRefs[i] == refIndex)
        {
            return true;
        }
    }
    return false;
}

YBI_INTEGRATOR_HD bool LightHasShadowExcludes(const LaunchParams &params, int lightIndex)
{
    const PackedLight *lights = GetPackedLights(params);
    return lights && lightIndex >= 0 && lightIndex < params.lightCount &&
           lights[lightIndex].shadowExcludeCount > 0u;
}

YBI_INTEGRATOR_HD bool PickDirectLightUniform(const LaunchParams &params,
                                              unsigned int &rngState,
                                              int *outLightIndex,
                                              float *outLightPickPdf)
{
    if (!outLightIndex || !outLightPickPdf || params.lights == 0ull || params.lightCount <= 0)
    {
        return false;
    }

    const int directCount = CountDirectLights(params);
    if (directCount <= 0)
    {
        return false;
    }

    const int target = ClampInt(int(Random01(rngState) * float(directCount)), 0, directCount - 1);
    const PackedLight *lights = GetPackedLights(params);
    int directIndex = 0;
    for (int i = 0; i < params.lightCount; ++i)
    {
        if (lights[i].type == static_cast<unsigned int>(LightType::Dome))
        {
            continue;
        }
        if (directIndex == target)
        {
            *outLightIndex = i;
            *outLightPickPdf = 1.0f / float(directCount);
            return true;
        }
        directIndex++;
    }
    return false;
}

YBI_INTEGRATOR_HD float SolidAnglePdfFromArea(float pickPdf,
                                              float area,
                                              float distanceSquared,
                                              float cosLight)
{
    if (pickPdf <= 0.0f || area <= 1.0e-8f || distanceSquared <= 1.0e-8f || cosLight <= 1.0e-8f)
    {
        return 0.0f;
    }
    return pickPdf * distanceSquared / (area * cosLight);
}

YBI_INTEGRATOR_HD bool FinalizeAreaLightSample(int lightIndex,
                                               float pickPdf,
                                               const PackedLight &light,
                                               const Vec3 &surfacePoint,
                                               const Vec3 &lightPoint,
                                               const Vec3 &lightNormal,
                                               DirectLightSample *outSample)
{
    if (!outSample)
    {
        return false;
    }

    const Vec3 toLight = Sub(lightPoint, surfacePoint);
    const float distanceSquared = Dot(toLight, toLight);
    if (distanceSquared <= 1.0e-8f)
    {
        return false;
    }

    const float distance = sqrtf(distanceSquared);
    const Vec3 wi = Mul(toLight, 1.0f / distance);
    const float cosLight = Dot(lightNormal, Mul(wi, -1.0f));
    if (cosLight <= 1.0e-8f)
    {
        return false;
    }

    outSample->lightIndex = lightIndex;
    outSample->wi = wi;
    outSample->radiance = LightEmission(light);
    outSample->lightPoint = lightPoint;
    outSample->lightNormal = lightNormal;
    outSample->distance = distance;
    outSample->pdf = SolidAnglePdfFromArea(pickPdf, light.areaScale, distanceSquared, cosLight);
    outSample->isDeltaLight = false;
    return outSample->pdf > 0.0f;
}

YBI_INTEGRATOR_HD bool SampleRectLight(int lightIndex,
                                       float pickPdf,
                                       const PackedLight &light,
                                       const Vec3 &surfacePoint,
                                       unsigned int &rngState,
                                       DirectLightSample *outSample)
{
    Vec3 tangent = {};
    Vec3 bitangent = {};
    Vec3 normal = {};
    GetLightBasis(light, &tangent, &bitangent, &normal);

    const float x = (Random01(rngState) - 0.5f) * light.width;
    const float y = (Random01(rngState) - 0.5f) * light.height;
    const Vec3 center = ToVec3(light.position);
    const Vec3 lightPoint = Add(Add(center, Mul(tangent, x)), Mul(bitangent, y));
    return FinalizeAreaLightSample(lightIndex, pickPdf, light, surfacePoint, lightPoint, normal, outSample);
}

YBI_INTEGRATOR_HD bool SampleDiskLight(int lightIndex,
                                       float pickPdf,
                                       const PackedLight &light,
                                       const Vec3 &surfacePoint,
                                       unsigned int &rngState,
                                       DirectLightSample *outSample)
{
    Vec3 tangent = {};
    Vec3 bitangent = {};
    Vec3 normal = {};
    GetLightBasis(light, &tangent, &bitangent, &normal);

    const float r = light.radius * SafeSqrt(Random01(rngState));
    const float phi = 6.28318530718f * Random01(rngState);
    const Vec3 localPoint = Add(Mul(tangent, r * cosf(phi)), Mul(bitangent, r * sinf(phi)));
    const Vec3 lightPoint = Add(ToVec3(light.position), localPoint);
    return FinalizeAreaLightSample(lightIndex, pickPdf, light, surfacePoint, lightPoint, normal, outSample);
}

YBI_INTEGRATOR_HD Vec3 SampleUniformSphereDirection(float u1, float u2)
{
    const float z = 1.0f - 2.0f * u1;
    const float r = SafeSqrt(1.0f - z * z);
    const float phi = 6.28318530718f * u2;
    return MakeVec3(r * cosf(phi), r * sinf(phi), z);
}

YBI_INTEGRATOR_HD bool SampleSphereLight(int lightIndex,
                                         float pickPdf,
                                         const PackedLight &light,
                                         const Vec3 &surfacePoint,
                                         unsigned int &rngState,
                                         DirectLightSample *outSample)
{
    const Vec3 normal = SampleUniformSphereDirection(Random01(rngState), Random01(rngState));
    const Vec3 lightPoint = Add(ToVec3(light.position), Mul(normal, light.radius));
    return FinalizeAreaLightSample(lightIndex, pickPdf, light, surfacePoint, lightPoint, normal, outSample);
}

YBI_INTEGRATOR_HD bool SampleCylinderLight(int lightIndex,
                                           float pickPdf,
                                           const PackedLight &light,
                                           const Vec3 &surfacePoint,
                                           unsigned int &rngState,
                                           DirectLightSample *outSample)
{
    Vec3 tangent = {};
    Vec3 bitangent = {};
    Vec3 axis = {};
    GetLightBasis(light, &tangent, &bitangent, &axis);

    const float phi = 6.28318530718f * Random01(rngState);
    const float z = (Random01(rngState) - 0.5f) * light.length;
    const Vec3 radial = Add(Mul(tangent, cosf(phi)), Mul(bitangent, sinf(phi)));
    const Vec3 center = ToVec3(light.position);
    const Vec3 lightPoint = Add(Add(center, Mul(radial, light.radius)), Mul(axis, z));
    return FinalizeAreaLightSample(lightIndex, pickPdf, light, surfacePoint, lightPoint, radial, outSample);
}

YBI_INTEGRATOR_HD bool SampleDirectLight(const LaunchParams &params,
                                         const Vec3 &surfacePoint,
                                         unsigned int &rngState,
                                         DirectLightSample *outSample)
{
    if (!outSample || params.lights == 0ull || params.lightCount <= 0)
    {
        return false;
    }

    int lightIndex = -1;
    float lightPickPdf = 0.0f;
    if (!PickDirectLightUniform(params, rngState, &lightIndex, &lightPickPdf))
    {
        return false;
    }

    const PackedLight *lights = GetPackedLights(params);
    const PackedLight &light = lights[lightIndex];
    switch (light.type)
    {
        case static_cast<unsigned int>(LightType::Distant):
            outSample->lightIndex = lightIndex;
            outSample->wi = Normalize(Mul(ToVec3(light.direction), -1.0f));
            outSample->radiance = LightEmission(light);
            outSample->distance = 1.0e20f;
            outSample->pdf = lightPickPdf;
            outSample->isDeltaLight = true;
            return true;
        case static_cast<unsigned int>(LightType::Rect):
            return SampleRectLight(lightIndex, lightPickPdf, light, surfacePoint, rngState, outSample);
        case static_cast<unsigned int>(LightType::Disk):
            return SampleDiskLight(lightIndex, lightPickPdf, light, surfacePoint, rngState, outSample);
        case static_cast<unsigned int>(LightType::Sphere):
            return SampleSphereLight(lightIndex, lightPickPdf, light, surfacePoint, rngState, outSample);
        case static_cast<unsigned int>(LightType::Cylinder):
            return SampleCylinderLight(lightIndex, lightPickPdf, light, surfacePoint, rngState, outSample);
        default:
            return false;
    }
}

#include "render/integrator_light_trace.inl"

} // namespace integrator
} // namespace render
} // namespace ybi
