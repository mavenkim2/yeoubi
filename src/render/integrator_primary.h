#pragma once

#include "render/integrator_texture.h"
#include "render/shading_core.h"

namespace ybi
{
namespace render
{
namespace integrator
{

template <typename State>
YBI_INTEGRATOR_HD uint32_t IntegratorPrimaryDiffuse(State &state, const HitInfo &hit)
{
    if (!hit.hasBarycentrics)
    {
        const uint32_t packed = ybi::render::PackRGB8(0.7f, 0.7f, 0.7f);
        return packed;
    }

    const LaunchParams &params = state.Params();
    const LaunchParams::InstanceGeomRef *refs =
        reinterpret_cast<const LaunchParams::InstanceGeomRef *>(params.instanceGeomRefs);
    if (!refs || hit.instanceId < 0 || hit.instanceId >= params.instanceGeomRefCount)
    {
        const uint32_t packed = ybi::render::PackRGB8(0.7f, 0.7f, 0.7f);
        return packed;
    }

    const LaunchParams::InstanceGeomRef geomRef = refs[hit.instanceId];
    Vec3 color = MakeVec3(0.7f, 0.7f, 0.7f);
    const bool sampled = TrySampleMaterialTexture(state, geomRef, hit, color);

    float outR = 0.0f;
    float outG = 0.0f;
    float outB = 0.0f;
    ybi::render::EvalPrimaryDiffuseCore(sampled, color.x, color.y, color.z, outR, outG, outB);
    return ybi::render::PackRGB8(outR, outG, outB);
}

template <typename State>
YBI_INTEGRATOR_HD void IntegratorFeedbackOnly(State &state, const HitInfo &hit)
{
    if (!hit.hasBarycentrics)
    {
        return;
    }

    const LaunchParams &params = state.Params();
    const LaunchParams::InstanceGeomRef *refs =
        reinterpret_cast<const LaunchParams::InstanceGeomRef *>(params.instanceGeomRefs);
    if (!refs || hit.instanceId < 0 || hit.instanceId >= params.instanceGeomRefCount)
    {
        return;
    }

    const LaunchParams::InstanceGeomRef geomRef = refs[hit.instanceId];
    TryWriteFeedbackOnly(state, geomRef, hit);
}

} // namespace integrator
} // namespace render
} // namespace ybi
