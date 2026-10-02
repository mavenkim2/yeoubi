#pragma once

#include "device/cuda_compat.h"
#include "kernel/bsdf/openpbr/ggx.h"

namespace ybi
{

YBI_GPU bool SampleThinWallTranslucentLobe(const OpenPbrBsdf &bsdf,
                                           const Vec2f &u,
                                           const Vec3f &viewDirection,
                                           Vec3f &lightDirection,
                                           Vec3f &weight)
{
    if (viewDirection.z <= 0.f || bsdf.thinWallTransmissionColor.MaxComponent() == 0.f)
        return false;

    bool success = SampleGGX(
        viewDirection,
        bsdf.specularRoughness,
        u,
        bsdf.specularReflectionIor,
        lightDirection,
        weight);
    if (!success)
        return false;

    // The sheet transmittance already includes Fresnel. Replace the reflection
    // weight with the constant transmission coefficient used by the reference.
    weight = bsdf.thinWallTransmissionColor *
             EvalAnisoSmithG1(lightDirection, bsdf.specularRoughness);
    lightDirection.z *= -1;
    return true;
}

YBI_GPU Vec3f EvaluateThinWallTranslucentLobe(const OpenPbrBsdf &bsdf,
                                              const Vec3f &viewDirection,
                                              const Vec3f &lightDirection)
{
    if (viewDirection.z <= 0.f || lightDirection.z >= 0.f ||
        bsdf.thinWallTransmissionColor.MaxComponent() == 0.f)
    {
        return Vec3f(0.f);
    }

    // In our face-forward local frame, the reference's reflected view and
    // flipped basis are equivalent to flipping the light direction's z.
    Vec3f reflectedLight(lightDirection.x, lightDirection.y, -lightDirection.z);
    float pdf = PDFGGX(viewDirection, reflectedLight, bsdf.specularRoughness);
    float shadowing = EvalAnisoSmithG1(reflectedLight, bsdf.specularRoughness);
    // BSDF * abs(cos(theta_l)) = color * D * G1(v) * G1(l) / (4 * v.z).
    return bsdf.thinWallTransmissionColor * (pdf * shadowing);
}

YBI_GPU float PDFThinWallTranslucentLobe(const OpenPbrBsdf &bsdf,
                                         const Vec3f &viewDirection,
                                         const Vec3f &lightDirection)
{
    if (viewDirection.z <= 0.f || lightDirection.z >= 0.f)
        return 0.f;

    Vec3f reflectedLight(lightDirection.x, lightDirection.y, -lightDirection.z);
    return PDFGGX(viewDirection, reflectedLight, bsdf.specularRoughness);
}

} // namespace ybi
