#pragma once

#include "render/bsdf/openpbr/ggx.h"
#include "render/bsdf/openpbr/glossy_diffuse.h"
#include "render/bsdf/openpbr/init.h"

namespace ybi
{

YBI_GPU bool SampleOpenPBR(const OpenPbrBsdf &bsdf,
                           const Vec3f &viewDirection,
                           const Vec3f &throughput,
                           Vec2f u,
                           Vec3f &lightDirection,
                           Vec3f &weight)
{
#if 0
    // Calculate specular lobe weight
    float Fr = FrDielectric(viewDirection.z, bsdf.specularIor);

    float specularLobeWeight = (Fr * throughput).MaxComponent();

    float diffuseLobeWeight = 0.f;

    float rand = u.x;

    float totalWeight = specularLobeWeight + diffuseLobeWeight;
    float specularLobeProbability = specularLobeWeight / totalWeight;
    float diffuseLobeProbability = 1.f - specularLobeProbability;
#endif

    bool success = true;
    int lobeIndex = 0;
    if (u.x < 0.5f) // diffuseLobeWeight)
    {
        lobeIndex = 0;
        u.x = u.x / 0.5f;
        success = SampleDiffuseLobe(viewDirection, bsdf.baseColor, u, lightDirection, weight);
    }
    else
    {
        lobeIndex = 1;
        u.x = (u.x - 0.5f) / (1 - 0.5f);
        success = SampleGlossLobe(viewDirection,
                                  bsdf.specularRoughness,
                                  u,
                                  bsdf.specularIor,
                                  bsdf.specularColor,
                                  lightDirection,
                                  weight);
    }

    if (!success)
        return false;

    // Calculate total pdf
    {
        float diffusePdf = OrenNayarPDF(lightDirection);
        float specularPdf = PDFGGX(viewDirection, lightDirection, bsdf.specularRoughness);

        float pdf = diffusePdf * 0.5f + specularPdf * 0.5f;

        if (pdf <= 0.f)
            return false;

        weight = Vec3f(0.f);
        // TODO: not hardcoded 0. also the diffuse lobe should be attenuated by gloss
        // reflection
        weight += EvaluateDiffuseLobe(viewDirection, lightDirection, bsdf.baseColor, 0.f);
        weight += EvaluateGlossLobe(viewDirection,
                                    lightDirection,
                                    bsdf.specularColor,
                                    bsdf.specularIor,
                                    bsdf.specularRoughness);

        weight /= pdf;
    }
    return true;

    // 1. select the lobe. seems it's just maxcomponent(fresnel * throughput)
    // 2. get the tables for diffuse to weight the light it receives
    // 3. profit
}

} // namespace ybi
