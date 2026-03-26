#pragma once

#include "render/integrator_common.h"
#include "render/launch_params.h"
#include "render/sampling/cdf_lookup.h"
#include "scene/material_light.h"

#include <cassert>

namespace ybi
{
namespace render
{
namespace integrator
{

struct DirectLightSample
{
    int lightIndex = -1;
    Vec3 wi = Vec3(0.0f, 0.0f, 1.0f);
    Vec3 radiance = Vec3(0.0f, 0.0f, 0.0f);
    Vec3 lightPoint = Vec3(0.0f, 0.0f, 0.0f);
    Vec3 lightNormal = Vec3(0.0f, 0.0f, 1.0f);
    float distance = 0.0f;
    float pdf = 0.0f;
    bool isDeltaLight = false;
};

struct LightRayHit
{
    int lightIndex = -1;
    Vec3 radiance = Vec3(0.0f, 0.0f, 0.0f);
    float distance = 0.0f;
    float pdf = 0.0f;
    bool isDeltaLight = false;
};

YBI_DEVICE float MaxF(float a, float b)
{
    return a > b ? a : b;
}

YBI_DEVICE float MinF(float a, float b)
{
    return a < b ? a : b;
}

YBI_DEVICE float SafeSqrtF(float x)
{
    return sqrtf(MaxF(x, 0.0f));
}

YBI_DEVICE Vec3 ToVec3(const ybi::Vec3 &value)
{
    return Vec3(value.x, value.y, value.z);
}

YBI_DEVICE const PackedLight *GetPackedLights(const LaunchParams &params);
YBI_DEVICE const int *GetLightShadowExcludeRefs(const LaunchParams &params);
YBI_DEVICE float DirectLightPickPdf(const LaunchParams &params, int lightIndex);
YBI_DEVICE bool
IsRefExcludedFromLightShadow(const LaunchParams &params, int lightIndex, int refIndex);
YBI_DEVICE bool LightHasShadowExcludes(const LaunchParams &params, int lightIndex);

YBI_DEVICE Vec3 LightEmission(const PackedLight &light)
{
    return ToVec3(light.color) * light.emissionScale;
}

YBI_DEVICE float LightSelectionWeight(const PackedLight &light)
{
    return light.selectionWeight > 0.0f ? light.selectionWeight : 0.0f;
}

YBI_DEVICE bool IsFiniteAreaLight(const PackedLight &light)
{
    return light.type == static_cast<unsigned int>(LightType::Rect) ||
           light.type == static_cast<unsigned int>(LightType::Disk) ||
           light.type == static_cast<unsigned int>(LightType::Sphere) ||
           light.type == static_cast<unsigned int>(LightType::Cylinder);
}

YBI_DEVICE Vec3 GetLightPosition(const PackedLight &light)
{
    return ybi::TransformPointAffine(light.worldFromLocal, Vec3(0.0f));
}

YBI_DEVICE Vec3 GetLightAxisXRaw(const PackedLight &light)
{
    return ybi::TransformVectorAffine(light.worldFromLocal, Vec3(1.0f, 0.0f, 0.0f));
}

YBI_DEVICE Vec3 GetLightAxisYRaw(const PackedLight &light)
{
    return ybi::TransformVectorAffine(light.worldFromLocal, Vec3(0.0f, 1.0f, 0.0f));
}

YBI_DEVICE Vec3 GetLightAxisZRaw(const PackedLight &light)
{
    return ybi::TransformVectorAffine(light.worldFromLocal, Vec3(0.0f, 0.0f, 1.0f));
}

YBI_DEVICE Vec3 SafeNormalizeAxis(const Vec3 &axis, const Vec3 &fallback)
{
    return LengthSquared(axis) > 1.0e-20f ? Normalize(axis) : fallback;
}

YBI_DEVICE void GetLightFrame(const PackedLight &light, Vec3 *outX, Vec3 *outY, Vec3 *outZ)
{
    const Vec3 z = SafeNormalizeAxis(GetLightAxisZRaw(light), Vec3(0.0f, 0.0f, 1.0f));
    Vec3 x = GetLightAxisXRaw(light);
    if (LengthSquared(x) <= 1.0e-20f)
    {
        x = Vec3(1.0f, 0.0f, 0.0f);
    }
    x = x - z * Dot(x, z);
    Vec3 y = Vec3(0.0f);
    if (LengthSquared(x) <= 1.0e-20f)
    {
        BuildOrthonormalBasis(z, x, y);
    }
    else
    {
        x = Normalize(x);
        y = Normalize(Cross(z, x));
        x = Normalize(Cross(y, z));
    }

    if (outX)
    {
        *outX = x;
    }
    if (outY)
    {
        *outY = y;
    }
    if (outZ)
    {
        *outZ = z;
    }
}

YBI_DEVICE Vec3 WorldDirectionToLightLocal(const PackedLight &light, const Vec3 &worldDirection)
{
    Vec3 x = {};
    Vec3 y = {};
    Vec3 z = {};
    GetLightFrame(light, &x, &y, &z);
    const Vec3 d = Normalize(worldDirection);
    return Vec3(Dot(d, x), Dot(d, y), Dot(d, z));
}

YBI_DEVICE int FindDomeLightIndex(const LaunchParams &params)
{
    const PackedLight *lights = GetPackedLights(params);
    if (!lights)
    {
        return -1;
    }

    for (int i = 0; i < params.lightCount; ++i)
    {
        if (lights[i].type == static_cast<unsigned int>(LightType::Dome))
        {
            return i;
        }
    }
    return -1;
}

YBI_DEVICE const char *LightTypeName(unsigned int type)
{
    switch (type)
    {
        case static_cast<unsigned int>(LightType::Dome):
            return "dome";
        case static_cast<unsigned int>(LightType::Distant):
            return "distant";
        case static_cast<unsigned int>(LightType::Rect):
            return "rect";
        case static_cast<unsigned int>(LightType::Disk):
            return "disk";
        case static_cast<unsigned int>(LightType::Sphere):
            return "sphere";
        case static_cast<unsigned int>(LightType::Cylinder):
            return "cylinder";
        default:
            return "unknown";
    }
}

YBI_DEVICE float DomeDirectionPdf()
{
    return 0.25f * ybi::kInvPi;
}

YBI_DEVICE void DirectionToLatLongUv(const Vec3 &direction, float *outU, float *outV);

YBI_DEVICE bool HasDomeCdf(const LaunchParams &params)
{
    return params.domeConditionalCdf != 0ull && params.domeMarginalCdf != 0ull &&
           params.domeCdfWidth > 0 && params.domeCdfHeight > 0;
}

YBI_DEVICE Vec3 LatLongUvToDirection(float u, float v)
{
    const float theta = ybi::kPi * Clamp(v, 0.0f, 1.0f);
    const float lon = (0.5f - u) * ybi::kTwoPi;
    const float sinTheta = sinf(theta);
    return Vec3(sinTheta * sinf(lon), cosf(theta), sinTheta * cosf(lon));
}

YBI_DEVICE float DomeTexelSolidAngle(int width, int height, int y)
{
    if (width <= 0 || height <= 0)
    {
        return 0.0f;
    }

    const int row = ClampInt(y, 0, height - 1);
    const float phiStep = ybi::kTwoPi / float(width);
    const float theta0 = ybi::kPi * (float(row) / float(height));
    const float theta1 = ybi::kPi * (float(row + 1) / float(height));
    return MaxF(phiStep * MaxF(cosf(theta0) - cosf(theta1), 0.0f), 0.0f);
}

YBI_DEVICE float DomeDirectionPdf(const LaunchParams &params,
                                  const PackedLight &domeLight,
                                  const Vec3 &worldDirection)
{
    if (!HasDomeCdf(params))
    {
        return DomeDirectionPdf();
    }

    const float *conditional = sampling::CdfPtr(params.domeConditionalCdf);
    const float *marginal = sampling::CdfPtr(params.domeMarginalCdf);
    if (!conditional || !marginal)
    {
        return DomeDirectionPdf();
    }

    float u = 0.5f;
    float v = 0.5f;
    const Vec3 localDirection = WorldDirectionToLightLocal(domeLight, worldDirection);
    DirectionToLatLongUv(localDirection, &u, &v);
    const int x = ClampInt(int(u * float(params.domeCdfWidth)), 0, params.domeCdfWidth - 1);
    const int y = ClampInt(int(v * float(params.domeCdfHeight)), 0, params.domeCdfHeight - 1);
    const float pmf = sampling::EvaluateCdf2DPmf(
        conditional, marginal, params.domeCdfWidth, params.domeCdfHeight, x, y);
    const float solidAngle = DomeTexelSolidAngle(params.domeCdfWidth, params.domeCdfHeight, y);
    return solidAngle > 0.0f ? pmf / solidAngle : 0.0f;
}

YBI_DEVICE bool IsDistantLight(const PackedLight &light)
{
    return light.type == static_cast<unsigned int>(LightType::Distant);
}

YBI_DEVICE bool IsFiniteAngleDistantLight(const PackedLight &light)
{
    return IsDistantLight(light) && light.solidAngle > 1.0e-8f;
}

YBI_DEVICE Vec3 LightLocalDirectionToWorld(const PackedLight &light, const Vec3 &localDirection)
{
    Vec3 x = {};
    Vec3 y = {};
    Vec3 z = {};
    GetLightFrame(light, &x, &y, &z);
    return Normalize(x * localDirection.x + y * localDirection.y + z * localDirection.z);
}

YBI_DEVICE Vec3 GetDistantLightDirection(const PackedLight &light)
{
    Vec3 axisX = {};
    Vec3 axisY = {};
    Vec3 axisZ = {};
    GetLightFrame(light, &axisX, &axisY, &axisZ);
    return axisZ;
}

YBI_DEVICE float DistantLightDirectionalPdf(const PackedLight &light)
{
    if (!IsFiniteAngleDistantLight(light))
    {
        return 0.0f;
    }
    return 1.0f / MaxF(light.solidAngle, 1.0e-8f);
}

YBI_DEVICE bool DirectionInsideDistantLightCap(const PackedLight &light, const Vec3 &worldDirection)
{
    if (!IsFiniteAngleDistantLight(light))
    {
        return false;
    }
    const Vec3 d = Normalize(worldDirection);
    return Dot(d, GetDistantLightDirection(light)) >= light.cosThetaMax;
}

YBI_DEVICE Vec3 EvaluateDistantLightRadiance(const PackedLight &light, const Vec3 &worldDirection)
{
    return DirectionInsideDistantLightCap(light, worldDirection) ? LightEmission(light)
                                                                 : Vec3(0.0f, 0.0f, 0.0f);
}

YBI_DEVICE Vec3 SampleUniformSphericalCapDirection(float cosThetaMax, float u1, float u2)
{
    const float z = Lerp(cosThetaMax, 1.0f, u1);
    const float r = SafeSqrtF(1.0f - z * z);
    const float phi = ybi::kTwoPi * u2;
    return Vec3(r * cosf(phi), r * sinf(phi), z);
}

YBI_DEVICE Vec3 SampleDistantLightDirection(const PackedLight &light, unsigned int &rngState)
{
    const Vec3 localDirection =
        SampleUniformSphericalCapDirection(light.cosThetaMax, Random01(rngState), Random01(rngState));
    return LightLocalDirectionToWorld(light, localDirection);
}

YBI_DEVICE void DirectionToLatLongUv(const Vec3 &direction, float *outU, float *outV)
{
    const Vec3 d = Normalize(direction);
    const float lon = atan2f(d.x, d.z);
    const float clampedY = MaxF(-1.0f, MinF(1.0f, d.y));
    if (outU)
    {
        const float u = 0.5f - lon * (0.5f * ybi::kInvPi);
        *outU = u - floorf(u);
    }
    if (outV)
    {
        *outV = acosf(clampedY) * ybi::kInvPi;
    }
}

template <typename State>
YBI_DEVICE Vec3 EvaluateDomeLightRadiance(State &state, const Vec3 &direction)
{
    const LaunchParams &params = state.Params();
    const int domeLightIndex = FindDomeLightIndex(params);
    const PackedLight *lights = GetPackedLights(params);
    const LaunchParams::MaterialTextureRef &textureRef = params.domeTextureRef;
    if (domeLightIndex < 0 || !lights)
    {
        return Vec3(0.0f, 0.0f, 0.0f);
    }
    const Vec3 domeEmission = LightEmission(lights[domeLightIndex]);
    if (textureRef.textureObject == 0ull || textureRef.valid == 0 || textureRef.width <= 0 ||
        textureRef.height <= 0)
    {
        return domeEmission;
    }

    float u = 0.5f;
    float v = 0.5f;
    Vec3 localDirection = direction;
    localDirection = WorldDirectionToLightLocal(lights[domeLightIndex], direction);
    DirectionToLatLongUv(localDirection, &u, &v);
    Vec4 sample = {};
    bool success = state.SampleTexture2D(textureRef, u, v, sample);
    assert(success);
    return Vec3(sample.x, sample.y, sample.z) * lights[domeLightIndex].emissionScale;
}

template <typename State>
YBI_DEVICE Vec3 EvaluateBackgroundRadiance(State &state, const Vec3 &direction)
{
    const LaunchParams &params = state.Params();
    return FindDomeLightIndex(params) >= 0 ? EvaluateDomeLightRadiance(state, direction)
                                           : SkyColor(direction);
}

template <typename State>
YBI_DEVICE Vec3 EvaluateMissInfiniteLightRadiance(State &state,
                                                  const Vec3 &direction,
                                                  bool hasBsdfContext,
                                                  bool skipNeeMis,
                                                  float bsdfPdf)
{
    const LaunchParams &params = state.Params();
    const PackedLight *lights = GetPackedLights(params);
    Vec3 radiance = EvaluateBackgroundRadiance(state, direction);

    const int domeLightIndex = FindDomeLightIndex(params);
    if (domeLightIndex >= 0 && hasBsdfContext && !skipNeeMis)
    {
        const float domePdf = DomeDirectionPdf(params, lights[domeLightIndex], direction);
        const float domeMis = bsdfPdf / MaxF(bsdfPdf + domePdf, 1.0e-6f);
        radiance *= domeMis;
    }

    if (!lights)
    {
        return radiance;
    }

    for (int i = 0; i < params.lightCount; ++i)
    {
        const PackedLight &light = lights[i];
        if (!IsFiniteAngleDistantLight(light))
        {
            continue;
        }

        const Vec3 lightRadiance = EvaluateDistantLightRadiance(light, direction);
        if (MaxComponent(lightRadiance) <= 0.0f)
        {
            continue;
        }

        float misWeight = 1.0f;
        if (hasBsdfContext && !skipNeeMis)
        {
            const float lightPdf = DirectLightPickPdf(params, i) * DistantLightDirectionalPdf(light);
            misWeight = bsdfPdf / MaxF(bsdfPdf + lightPdf, 1.0e-6f);
        }
        radiance += lightRadiance * misWeight;
    }

    return radiance;
}

YBI_DEVICE int CountDirectLights(const LaunchParams &params)
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

YBI_DEVICE float DirectLightPickPdf(const LaunchParams &params)
{
    const int count = CountDirectLights(params);
    return count > 0 ? 1.0f / float(count) : 0.0f;
}

YBI_DEVICE float TotalDirectLightSelectionWeight(const LaunchParams &params)
{
    const PackedLight *lights = GetPackedLights(params);
    if (!lights)
    {
        return 0.0f;
    }

    float total = 0.0f;
    for (int i = 0; i < params.lightCount; ++i)
    {
        if (lights[i].type == static_cast<unsigned int>(LightType::Dome))
        {
            continue;
        }
        total += LightSelectionWeight(lights[i]);
    }
    return total;
}

YBI_DEVICE float DirectLightPickPdf(const LaunchParams &params, int lightIndex)
{
    const PackedLight *lights = GetPackedLights(params);
    if (!lights || lightIndex < 0 || lightIndex >= params.lightCount ||
        lights[lightIndex].type == static_cast<unsigned int>(LightType::Dome))
    {
        return 0.0f;
    }

    const float totalWeight = TotalDirectLightSelectionWeight(params);
    if (totalWeight > 0.0f)
    {
        return LightSelectionWeight(lights[lightIndex]) / totalWeight;
    }
    return DirectLightPickPdf(params);
}

YBI_DEVICE const PackedLight *GetPackedLights(const LaunchParams &params)
{
    if (params.lights == 0ull || params.lightCount <= 0)
    {
        return nullptr;
    }
    return reinterpret_cast<const PackedLight *>(params.lights);
}

YBI_DEVICE const int *GetLightShadowExcludeRefs(const LaunchParams &params)
{
    if (params.lightShadowExcludeRefs == 0ull || params.lightShadowExcludeRefCount <= 0)
    {
        return nullptr;
    }
    return reinterpret_cast<const int *>(params.lightShadowExcludeRefs);
}

YBI_DEVICE bool
IsRefExcludedFromLightShadow(const LaunchParams &params, int lightIndex, int refIndex)
{
    const PackedLight *lights = GetPackedLights(params);
    const int *excludeRefs = GetLightShadowExcludeRefs(params);
    if (!lights || !excludeRefs || lightIndex < 0 || lightIndex >= params.lightCount ||
        refIndex < 0)
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

YBI_DEVICE bool LightHasShadowExcludes(const LaunchParams &params, int lightIndex)
{
    const PackedLight *lights = GetPackedLights(params);
    return lights && lightIndex >= 0 && lightIndex < params.lightCount &&
           lights[lightIndex].shadowExcludeCount > 0u;
}

YBI_DEVICE bool PickDirectLightWeighted(const LaunchParams &params,
                                        unsigned int &rngState,
                                        int *outLightIndex,
                                        float *outLightPickPdf)
{
    if (!outLightIndex || !outLightPickPdf || params.lights == 0ull || params.lightCount <= 0)
    {
        return false;
    }

    const PackedLight *lights = GetPackedLights(params);
    if (!lights)
    {
        return false;
    }

    const float totalWeight = TotalDirectLightSelectionWeight(params);
    if (totalWeight <= 0.0f)
    {
        const int directCount = CountDirectLights(params);
        if (directCount <= 0)
        {
            return false;
        }
        const int target =
            ClampInt(int(Random01(rngState) * float(directCount)), 0, directCount - 1);
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

    float targetWeight = Random01(rngState) * totalWeight;
    for (int i = 0; i < params.lightCount; ++i)
    {
        if (lights[i].type == static_cast<unsigned int>(LightType::Dome))
        {
            continue;
        }
        const float weight = LightSelectionWeight(lights[i]);
        if (weight <= 0.0f)
        {
            continue;
        }
        if (targetWeight <= weight)
        {
            *outLightIndex = i;
            *outLightPickPdf = weight / totalWeight;
            return true;
        }
        targetWeight -= weight;
    }
    for (int i = params.lightCount - 1; i >= 0; --i)
    {
        if (lights[i].type == static_cast<unsigned int>(LightType::Dome))
        {
            continue;
        }
        const float weight = LightSelectionWeight(lights[i]);
        if (weight > 0.0f)
        {
            *outLightIndex = i;
            *outLightPickPdf = weight / totalWeight;
            return true;
        }
    }
    return false;
}

YBI_DEVICE float
SolidAnglePdfFromArea(float pickPdf, float area, float distanceSquared, float cosLight)
{
    if (pickPdf <= 0.0f || area <= 1.0e-8f || distanceSquared <= 1.0e-8f || cosLight <= 1.0e-8f)
    {
        return 0.0f;
    }
    return pickPdf * distanceSquared / (area * cosLight);
}

YBI_DEVICE bool FinalizeAreaLightSample(int lightIndex,
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

    const Vec3 toLight = lightPoint - surfacePoint;
    const float distanceSquared = Dot(toLight, toLight);
    if (distanceSquared <= 1.0e-8f)
    {
        return false;
    }

    const float distance = sqrtf(distanceSquared);
    const Vec3 wi = toLight * (1.0f / distance);
    const float cosLight = Dot(lightNormal, -wi);
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

YBI_DEVICE bool SampleRectLight(int lightIndex,
                                float pickPdf,
                                const PackedLight &light,
                                const Vec3 &surfacePoint,
                                unsigned int &rngState,
                                DirectLightSample *outSample)
{
    Vec3 tangent = {};
    Vec3 bitangent = {};
    Vec3 normal = {};
    GetLightFrame(light, &tangent, &bitangent, &normal);
    normal = -normal;

    const float x = (Random01(rngState) - 0.5f) * light.width;
    const float y = (Random01(rngState) - 0.5f) * light.height;
    const Vec3 center = GetLightPosition(light);
    const Vec3 lightPoint = center + tangent * x + bitangent * y;
    return FinalizeAreaLightSample(
        lightIndex, pickPdf, light, surfacePoint, lightPoint, normal, outSample);
}

YBI_DEVICE bool SampleDiskLight(int lightIndex,
                                float pickPdf,
                                const PackedLight &light,
                                const Vec3 &surfacePoint,
                                unsigned int &rngState,
                                DirectLightSample *outSample)
{
    Vec3 tangent = {};
    Vec3 bitangent = {};
    Vec3 normal = {};
    GetLightFrame(light, &tangent, &bitangent, &normal);
    normal = -normal;

    const float r = light.radius * SafeSqrtF(Random01(rngState));
    const float phi = ybi::kTwoPi * Random01(rngState);
    const Vec3 localPoint = tangent * (r * cosf(phi)) + bitangent * (r * sinf(phi));
    const Vec3 lightPoint = GetLightPosition(light) + localPoint;
    return FinalizeAreaLightSample(
        lightIndex, pickPdf, light, surfacePoint, lightPoint, normal, outSample);
}

YBI_DEVICE Vec3 SampleUniformSphereDirection(float u1, float u2)
{
    const float z = 1.0f - 2.0f * u1;
    const float r = SafeSqrtF(1.0f - z * z);
    const float phi = ybi::kTwoPi * u2;
    return Vec3(r * cosf(phi), r * sinf(phi), z);
}

YBI_DEVICE bool SampleSphereLight(int lightIndex,
                                  float pickPdf,
                                  const PackedLight &light,
                                  const Vec3 &surfacePoint,
                                  unsigned int &rngState,
                                  DirectLightSample *outSample)
{
    const Vec3 normal = SampleUniformSphereDirection(Random01(rngState), Random01(rngState));
    const Vec3 lightPoint = GetLightPosition(light) + normal * light.radius;
    return FinalizeAreaLightSample(
        lightIndex, pickPdf, light, surfacePoint, lightPoint, normal, outSample);
}

YBI_DEVICE bool SampleCylinderLight(int lightIndex,
                                    float pickPdf,
                                    const PackedLight &light,
                                    const Vec3 &surfacePoint,
                                    unsigned int &rngState,
                                    DirectLightSample *outSample)
{
    Vec3 tangent = {};
    Vec3 bitangent = {};
    Vec3 axis = {};
    GetLightFrame(light, &tangent, &bitangent, &axis);

    const float phi = ybi::kTwoPi * Random01(rngState);
    const float z = (Random01(rngState) - 0.5f) * light.length;
    const Vec3 radial = tangent * cosf(phi) + bitangent * sinf(phi);
    const Vec3 center = GetLightPosition(light);
    const Vec3 lightPoint = center + radial * light.radius + axis * z;
    return FinalizeAreaLightSample(
        lightIndex, pickPdf, light, surfacePoint, lightPoint, radial, outSample);
}

YBI_DEVICE bool SampleDirectLight(const LaunchParams &params,
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
    if (!PickDirectLightWeighted(params, rngState, &lightIndex, &lightPickPdf))
    {
        return false;
    }

    const PackedLight *lights = GetPackedLights(params);
    const PackedLight &light = lights[lightIndex];
    switch (light.type)
    {
        case static_cast<unsigned int>(LightType::Distant):
        {
            outSample->lightIndex = lightIndex;
            outSample->radiance = LightEmission(light);
            outSample->distance = 1.0e20f;
            if (IsFiniteAngleDistantLight(light))
            {
                outSample->wi = SampleDistantLightDirection(light, rngState);
                outSample->pdf = lightPickPdf * DistantLightDirectionalPdf(light);
                outSample->isDeltaLight = false;
                return outSample->pdf > 0.0f;
            }

            outSample->wi = GetDistantLightDirection(light);
            outSample->pdf = lightPickPdf;
            outSample->isDeltaLight = true;
            return true;
        }
        case static_cast<unsigned int>(LightType::Rect):
            return SampleRectLight(
                lightIndex, lightPickPdf, light, surfacePoint, rngState, outSample);
        case static_cast<unsigned int>(LightType::Disk):
            return SampleDiskLight(
                lightIndex, lightPickPdf, light, surfacePoint, rngState, outSample);
        case static_cast<unsigned int>(LightType::Sphere):
            return SampleSphereLight(
                lightIndex, lightPickPdf, light, surfacePoint, rngState, outSample);
        case static_cast<unsigned int>(LightType::Cylinder):
            return SampleCylinderLight(
                lightIndex, lightPickPdf, light, surfacePoint, rngState, outSample);
        default:
            return false;
    }
}

template <typename State>
YBI_DEVICE bool SampleDomeLight(State &state, unsigned int &rngState, DirectLightSample *outSample)
{
    if (!outSample)
    {
        return false;
    }

    const LaunchParams &params = state.Params();
    const int domeLightIndex = FindDomeLightIndex(params);
    if (domeLightIndex < 0)
    {
        return false;
    }

    outSample->lightIndex = domeLightIndex;
    if (HasDomeCdf(params))
    {
        const float *conditional = sampling::CdfPtr(params.domeConditionalCdf);
        const float *marginal = sampling::CdfPtr(params.domeMarginalCdf);
        sampling::CdfSample2D cdfSample = {};
        if (conditional && marginal &&
            sampling::SampleCdf2D(conditional,
                                  marginal,
                                  params.domeCdfWidth,
                                  params.domeCdfHeight,
                                  Random01(rngState),
                                  Random01(rngState),
                                  &cdfSample))
        {
            const Vec3 localDirection = LatLongUvToDirection(cdfSample.u, cdfSample.v);
            outSample->wi = LightLocalDirectionToWorld(GetPackedLights(params)[domeLightIndex], localDirection);
            const float solidAngle =
                DomeTexelSolidAngle(params.domeCdfWidth, params.domeCdfHeight, cdfSample.y);
            outSample->pdf = solidAngle > 0.0f ? cdfSample.pmf / solidAngle : 0.0f;
        }
        else
        {
            outSample->wi = SampleUniformSphereDirection(Random01(rngState), Random01(rngState));
            outSample->pdf = DomeDirectionPdf();
        }
    }
    else
    {
        outSample->wi = SampleUniformSphereDirection(Random01(rngState), Random01(rngState));
        outSample->pdf = DomeDirectionPdf();
    }
    outSample->radiance = EvaluateDomeLightRadiance(state, outSample->wi);
    outSample->lightPoint = Vec3(0.0f, 0.0f, 0.0f);
    outSample->lightNormal = -outSample->wi;
    outSample->distance = 1.0e20f;
    outSample->isDeltaLight = false;
    return outSample->pdf > 0.0f && MaxComponent(outSample->radiance) > 0.0f;
}

#include "render/integrator_light_trace.inl"

} // namespace integrator
} // namespace render
} // namespace ybi
