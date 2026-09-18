#pragma once

#include "core/bit.h"
#include "core/math.h"
#include "render/bsdf/materialx/closure.h"
#include "render/bsdf/materialx/material_defines.h"
#include "render/bsdf/materialx/sheen_types.h"
#include "scene/device_scene.h"

namespace ybi
{

YBI_GPU void GetSheenParameters(const MaterialDescriptor &descriptor,
                                const BsdfClosure &closure,
                                SheenBSDF &bsdf)
{
    const int colorIndex =
        ExtractBits(descriptor.bits, MATERIAL_NUM_CLOSURE_BITS, SHEEN_TINT_OFFSET);
    const int roughnessIndex =
        ExtractBits(descriptor.bits, MATERIAL_NUM_CLOSURE_BITS, SHEEN_ROUGHNESS_OFFSET);
    const int probIndex = ExtractBits(descriptor.bits, 5, 4);

    bsdf.roughness = closure.values[roughnessIndex];
    // bsdf.probIndex = probIndex;
}

// NOTE: taken from materialx implementation
YBI_GPU float ImageWorksSheenDirectionalAlbedoAnalytic(float nDotV, float roughness)
{
    Vec2f r = Vec2f(13.67300, 1.0) + Vec2f(-68.78018, 61.57746) * nDotV +
              Vec2f(799.08825, 442.78211) * roughness +
              Vec2f(-905.00061, 2597.49308) * nDotV * roughness +
              Vec2f(60.28956, 121.81241) * Sqr(nDotV) +
              Vec2f(1086.96473, 3045.55075) * Sqr(roughness);
    return r.x / r.y;
}

YBI_GPU float SheenAlbedo(const SheenBSDF &bsdf, const Vec3f &viewDirection)
{
    return ImageWorksSheenDirectionalAlbedoAnalytic(viewDirection.z, bsdf.roughness);
}

YBI_GPU float SheenProbability(const SheenBSDF &sheen, const Vec3f &viewDirection)
{
    // ImageWorksSheenDirectionalAlbedoAnalytic(viewDirection.z, sheen.roughness);
    return 0.f;
}

// Conty kulla sheen
} // namespace ybi
