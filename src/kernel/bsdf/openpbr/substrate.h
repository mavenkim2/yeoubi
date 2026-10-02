#pragma once

#include "kernel/bsdf/openpbr/gloss.h"
#include "kernel/bsdf/openpbr/glossy_diffuse.h"
#include "kernel/bsdf/openpbr/translucent.h"

namespace ybi
{

YBI_GPU Vec3f EvaluateOpenPBRBaseSubstrate(const OpenPbrBsdf &bsdf,
                                           const Vec3f &viewDirection,
                                           const Vec3f &lightDirection)
{
    Vec3f weight(0.f);
    weight += EvaluateThinWallTranslucentLobe(bsdf, viewDirection, lightDirection);
    // TODO: not hardcoded 0. also the diffuse lobe should be attenuated by gloss
    // reflection

    weight += EvaluateDiffuseLobe(viewDirection, lightDirection, bsdf.baseColor, 0.f);
    weight += EvaluateGlossLobe(viewDirection,
                                lightDirection,
                                bsdf.specularColor,
                                bsdf.specularReflectionIor,
                                bsdf.specularRoughness);
    return weight;
}

YBI_GPU float PDFOpenPBRBaseSubstrate(const OpenPbrBsdf &bsdf,
                                      const Vec3f &viewDirection,
                                      const Vec3f &lightDirection)
{
    float diffusePdf = OrenNayarPDF(lightDirection);
    float specularPdf = PDFGGX(viewDirection, lightDirection, bsdf.specularRoughness);

    float thinWallTransmissionPdf =
        PDFThinWallTranslucentLobe(bsdf, viewDirection, lightDirection);
    const float *lobeWeights = bsdf.lobeWeights;
    float lobeWeightsSum =
        lobeWeights[THIN_WALL_TRANS] + lobeWeights[DIFFUSE] + lobeWeights[SPECULAR];
    if (lobeWeightsSum <= 0.f)
        return 0.f;

    float pdf = (thinWallTransmissionPdf * lobeWeights[THIN_WALL_TRANS] +
                 diffusePdf * lobeWeights[DIFFUSE] + specularPdf * lobeWeights[SPECULAR]) /
                lobeWeightsSum;

    pdf = Max(0.f, pdf);

    return pdf;
}

YBI_GPU bool SampleOpenPBRBaseSubstrate(const OpenPbrBsdf &bsdf,
                                        const Vec3f &viewDirection,
                                        const Vec3f &throughput,
                                        Vec2f u,
                                        Vec3f &lightDirection,
                                        Vec3f &weight)
{
    const float *lobeWeights = bsdf.lobeWeights;
    float lobeWeightsSum =
        lobeWeights[THIN_WALL_TRANS] + lobeWeights[DIFFUSE] + lobeWeights[SPECULAR];
    if (lobeWeightsSum <= 0.f)
        return false;
    float rand = u.x * lobeWeightsSum;
    float lobeWeight = 0.f;

    int selectedLobeIndex = 0;
    for (int lobeIndex = 0; lobeIndex < COUNT; lobeIndex++)
    {
        float selectedLobeWeight = lobeWeights[lobeIndex];
        if (rand < lobeWeight + selectedLobeWeight || lobeIndex == COUNT - 1)
        {
            selectedLobeIndex = lobeIndex;
            rand = (rand - lobeWeight) / selectedLobeWeight;
            rand = ClampRandomNumber(rand);
            break;
        }
        lobeWeight += lobeWeights[lobeIndex];
    }

    u = Vec2f(rand, u.y);
    bool success = true;
    if (selectedLobeIndex == THIN_WALL_TRANS)
    {
        success = SampleThinWallTranslucentLobe(bsdf, u, viewDirection, lightDirection, weight);
    }
    else if (selectedLobeIndex == DIFFUSE)
    {
        success = SampleDiffuseLobe(viewDirection, bsdf.baseColor, u, lightDirection, weight);
    }
    else if (selectedLobeIndex == SPECULAR)
    {
        success = SampleGlossLobe(viewDirection,
                                  bsdf.specularRoughness,
                                  u,
                                  bsdf.specularReflectionIor,
                                  bsdf.specularColor,
                                  lightDirection,
                                  weight);
    }

    if (!success)
        return false;

    float pdf = PDFOpenPBRBaseSubstrate(bsdf, viewDirection, lightDirection);
    if (pdf <= 0.f)
        return false;

    weight = EvaluateOpenPBRBaseSubstrate(bsdf, viewDirection, lightDirection) / pdf;
    return true;
}

} // namespace ybi
