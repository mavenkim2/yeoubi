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
    Vec3 wi = MakeVec3(0.0f, 0.0f, 1.0f);
    Vec3 radiance = MakeVec3(0.0f, 0.0f, 0.0f);
    float distance = 0.0f;
    float pdf = 0.0f;
};

YBI_INTEGRATOR_HD float MaxF(float a, float b)
{
    return a > b ? a : b;
}

YBI_INTEGRATOR_HD Vec3 ToVec3(const ybi::float3 &value)
{
    return MakeVec3(value.x, value.y, value.z);
}

YBI_INTEGRATOR_HD Vec3 EnvironmentRadiance(const LaunchParams &params, const Vec3 &direction)
{
    Vec3 env = MakeVec3(0.0f, 0.0f, 0.0f);
    bool hasDome = false;
    if (params.lights != 0ull && params.lightCount > 0)
    {
        const PackedLight *lights = reinterpret_cast<const PackedLight *>(params.lights);
        for (int i = 0; i < params.lightCount; ++i)
        {
            const PackedLight &light = lights[i];
            if (light.type != static_cast<unsigned int>(LightType::Dome))
            {
                continue;
            }
            env = Add(env, Mul(ToVec3(light.color), light.intensityScale));
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
    if (params.lights == 0ull || params.lightCount <= 0)
    {
        return 0;
    }

    int count = 0;
    const PackedLight *lights = reinterpret_cast<const PackedLight *>(params.lights);
    for (int i = 0; i < params.lightCount; ++i)
    {
        if (lights[i].type != static_cast<unsigned int>(LightType::Dome))
        {
            count++;
        }
    }
    return count;
}

YBI_INTEGRATOR_HD float DirectLightWeight(const PackedLight &light)
{
    if (light.type == static_cast<unsigned int>(LightType::Dome))
    {
        return 0.0f;
    }
    return light.selectionWeight > 0.0f ? light.selectionWeight : 0.0f;
}

YBI_INTEGRATOR_HD float TotalDirectLightWeight(const LaunchParams &params)
{
    if (params.lights == 0ull || params.lightCount <= 0)
    {
        return 0.0f;
    }

    float total = 0.0f;
    const PackedLight *lights = reinterpret_cast<const PackedLight *>(params.lights);
    for (int i = 0; i < params.lightCount; ++i)
    {
        total += DirectLightWeight(lights[i]);
    }
    return total;
}

YBI_INTEGRATOR_HD bool SampleDirectLight(const LaunchParams &params,
                                         const Vec3 &surfacePoint,
                                         unsigned int &rngState,
                                         DirectLightSample *outSample)
{
    if (!outSample)
    {
        return false;
    }

    const PackedLight *lights = reinterpret_cast<const PackedLight *>(params.lights);
    const float totalWeight = TotalDirectLightWeight(params);
    const int directCount = CountDirectLights(params);
    if (directCount <= 0)
    {
        return false;
    }

    const PackedLight *chosen = nullptr;
    float chosenPdf = 0.0f;
    if (totalWeight > 0.0f)
    {
        float target = Random01(rngState) * totalWeight;
        float accum = 0.0f;
        for (int i = 0; i < params.lightCount; ++i)
        {
            const float weight = DirectLightWeight(lights[i]);
            if (weight <= 0.0f)
            {
                continue;
            }
            accum += weight;
            if (target <= accum || i == params.lightCount - 1)
            {
                chosen = &lights[i];
                chosenPdf = weight / MaxF(totalWeight, 1.0e-6f);
                break;
            }
        }
    }
    else
    {
        const float randomIndex = Random01(rngState) * float(directCount);
        const int target = ClampInt(int(randomIndex), 0, directCount - 1);
        int directIndex = 0;
        for (int i = 0; i < params.lightCount; ++i)
        {
            if (lights[i].type == static_cast<unsigned int>(LightType::Dome))
            {
                continue;
            }
            if (directIndex == target)
            {
                chosen = &lights[i];
                chosenPdf = 1.0f / float(directCount);
                break;
            }
            directIndex++;
        }
    }
    if (!chosen)
    {
        return false;
    }

    outSample->pdf = chosenPdf;
    switch (chosen->type)
    {
        case static_cast<unsigned int>(LightType::Distant):
        {
            outSample->wi = Normalize(Mul(ToVec3(chosen->direction), -1.0f));
            outSample->distance = 1.0e20f;
            outSample->radiance = Mul(ToVec3(chosen->color), chosen->intensityScale);
            return true;
        }
        case static_cast<unsigned int>(LightType::Rect):
        case static_cast<unsigned int>(LightType::Disk):
        case static_cast<unsigned int>(LightType::Sphere):
        case static_cast<unsigned int>(LightType::Cylinder):
        {
            const Vec3 toLight = Sub(ToVec3(chosen->position), surfacePoint);
            const float dist2 = Dot(toLight, toLight);
            if (dist2 <= 1.0e-8f)
            {
                return false;
            }
            const float dist = sqrtf(dist2);
            const Vec3 wi = Mul(toLight, 1.0f / dist);
            float orientation = 1.0f;
            if ((chosen->flags & LIGHT_FLAG_ONE_SIDED) != 0u)
            {
                orientation = MaxF(0.0f, Dot(ToVec3(chosen->direction), Mul(wi, -1.0f)));
                if (orientation <= 0.0f)
                {
                    return false;
                }
            }

            outSample->wi = wi;
            outSample->distance = dist;
            outSample->radiance =
                Mul(ToVec3(chosen->color),
                    (chosen->intensityScale * chosen->areaScale * orientation) / dist2);
            return true;
        }
        default:
            return false;
    }
}

} // namespace integrator
} // namespace render
} // namespace ybi
