#pragma once

#include "core/math.h"
#include "core/vec2f.h"
#include "core/vec3f.h"
#include "render/bsdf/materialx/closure.h"
#include "render/bsdf/materialx/dielectric_bsdf.h"
#include "scene/device_scene.h"
#include <cstdint>

namespace ybi
{

enum ConductorParameterOffset : uint32_t
{
    ConductorWeight = 0,
    ConductorIor = 1,
    ConductorExtinction = 4,
    ConductorRoughness = 7,
    ConductorNormal = 9,
    ConductorParameterCount = 12,
};

struct ConductorBsdf
{
    float weight;
    Vec3f ior;
    Vec3f extinction;
    Vec2f alpha;
    Vec3f normal;
};

YBI_CPU_GPU bool GetConductorParameters(const MaterialDescriptor &descriptor,
                                        const BsdfClosure &closure,
                                        ConductorBsdf &bsdf)
{
    const uint32_t base = descriptor.parameterOffset;
    if (!closure.CanLoad(base, ConductorParameterCount))
        return false;
    return closure.LoadFloat(base + ConductorWeight, bsdf.weight) &&
           closure.LoadFloat3(base + ConductorIor, bsdf.ior) &&
           closure.LoadFloat3(base + ConductorExtinction, bsdf.extinction) &&
           closure.LoadFloat2(base + ConductorRoughness, bsdf.alpha) &&
           closure.LoadFloat3(base + ConductorNormal, bsdf.normal);
}

// A directional-albedo estimate for lobe importance, not a conductor BRDF.
YBI_CPU_GPU Vec3f ConductorDirectionalAlbedo(const ConductorBsdf &bsdf, float nDotV)
{
    if (bsdf.weight <= 0.f)
        return Vec3f(0.f);

    Vec3f f0;
    for (int channel = 0; channel < 3; ++channel)
    {
        const float n = bsdf.ior[channel];
        const float k2 = Sqr(bsdf.extinction[channel]);
        f0[channel] = (Sqr(n - 1.f) + k2) / (Sqr(n + 1.f) + k2);
    }

    const float alpha = Sqrt(Max(1e-8f, bsdf.alpha.x) * Max(1e-8f, bsdf.alpha.y));
    const Vec3f fAvg = f0 + (Vec3f(1.f) - f0) * (1.f / 21.f);
    const Vec3f comp = GGXEnergyCompensation(nDotV, alpha, fAvg);
    return Clamp(bsdf.weight * DielectricBsdfAlbedo(nDotV, alpha, f0, Vec3f(1.f)) * comp,
                 0.f,
                 1.f);
}

} // namespace ybi
