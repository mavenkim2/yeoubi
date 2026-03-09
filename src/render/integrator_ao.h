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
YBI_INTEGRATOR_HD uint32_t IntegratorAO(State &state, const HitInfo &hit)
{
    const LaunchParams &params = state.Params();
    const int spp = params.spp > 0 ? params.spp : 1;

    Vec3 normal = Vec3(-hit.rayDir.x, -hit.rayDir.y, -hit.rayDir.z);
    if (hit.hasGeomNormal)
    {
        normal = FaceForward(Normalize(hit.geomNormal), hit.rayDir);
    }
    Vec3 tangent = {};
    Vec3 bitangent = {};
    BuildOrthonormalBasis(normal, tangent, bitangent);

    int visible = 0;
    UInt2 launchIndex = state.LaunchIndex();
    unsigned int rngState = Hash32((launchIndex.x + 1u) * 73856093u ^
                                   (launchIndex.y + 1u) * 19349663u ^
                                   (static_cast<unsigned int>(params.currentSpp) + 1u) * 83492791u);

    for (int i = 0; i < spp; ++i)
    {
        const float u1 = Random01(rngState);
        const float u2 = Random01(rngState);
        const Vec3 local = SampleCosineHemisphere(u1, u2);
        Vec3 sampleDir = Vec3(tangent.x * local.x + bitangent.x * local.y + normal.x * local.z,
                                  tangent.y * local.x + bitangent.y * local.y + normal.y * local.z,
                                  tangent.z * local.x + bitangent.z * local.y + normal.z * local.z);
        sampleDir = Normalize(sampleDir);
        const Vec3 sampleOrigin = Vec3(hit.rayOrigin.x + hit.rayDir.x * hit.t +
                                               normal.x * params.aoBias,
                                           hit.rayOrigin.y + hit.rayDir.y * hit.t +
                                               normal.y * params.aoBias,
                                           hit.rayOrigin.z + hit.rayDir.z * hit.t +
                                               normal.z * params.aoBias);
        const bool occluded =
            state.TraceOcclusion(sampleOrigin, sampleDir, params.aoBias, params.aoMaxDistance);
        if (!occluded)
        {
            visible++;
        }
    }

    const float ao = ybi::render::EvalAOVisibility(visible, spp);
    return ybi::render::PackRGB8(ao, ao, ao);
}

} // namespace integrator
} // namespace render
} // namespace ybi
