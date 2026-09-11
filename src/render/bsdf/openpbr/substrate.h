#pragma once

#include "render/bsdf/openpbr/fuzz.h"
#include "render/bsdf/openpbr/ggx.h"
#include "render/bsdf/openpbr/glossy_diffuse.h"
#include "render/bsdf/openpbr/init.h"

namespace ybi
{

YBI_GPU bool SampleOpenPBRBaseSubstrate(const OpenPbrBsdf &bsdf,
                                        const Vec3f &viewDirection,
                                        const Vec3f &throughput,
                                        Vec2f u,
                                        Vec3f &lightDirection,
                                        Vec3f &weight);

YBI_GPU Vec3f EvaluateOpenPBRBaseSubstrate(const OpenPbrBsdf &bsdf,
                                           const Vec3f &viewDirection,
                                           const Vec3f &lightDirection);

YBI_GPU float PDFOpenPBRBaseSubstrate(const OpenPbrBsdf &bsdf,
                                      const Vec3f &viewDirection,
                                      const Vec3f &lightDirection);

YBI_GPU bool SampleOpenPBR(const OpenPbrBsdf &bsdf,
                           const Vec3f &viewDirection,
                           const Vec3f &throughput,
                           Vec2f u,
                           Vec3f &lightDirection,
                           Vec3f &weight)
{
    Vec3f ltcCoefficients = GetLTCCoefficients(bsdf.fuzzRoughness, viewDirection.z);
    float fuzzReflectance = ltcCoefficients[2];

    // Calculate sheen lobe weight
    float sheenWeight = bsdf.fuzzWeight * fuzzReflectance * bsdf.fuzzColor.MaxComponent();
    float sheenProb = sheenWeight / (1.f + sheenWeight);

#if 0
    // Calculate specular lobe weight float Fr = FrDielectric(viewDirection.z, bsdf.specularIor);
    float specularLobeWeight = (Fr * throughput).MaxComponent();

    float diffuseLobeWeight = 0.f;

    float rand = u.x;

    float totalWeight = specularLobeWeight + diffuseLobeWeight;
    float specularLobeProbability = specularLobeWeight / totalWeight;
    float diffuseLobeProbability = 1.f - specularLobeProbability;
#endif

    // Sample fuzz layer
    if (u.x < sheenProb)
    {
        u.x = u.x / sheenProb;

        bool success = SampleSheen(viewDirection, u, bsdf.fuzzRoughness, lightDirection);
        if (!success)
            return false;
    }
    else
    {
        u.x = (u.x - sheenProb) / (1 - sheenProb);

        bool success =
            SampleOpenPBRBaseSubstrate(bsdf, viewDirection, throughput, u, lightDirection, weight);

        if (!success)
            return false;
    }

    Vec3f substrateWeight = EvaluateOpenPBRBaseSubstrate(bsdf, viewDirection, lightDirection);
    float substratePdf = PDFOpenPBRBaseSubstrate(bsdf, viewDirection, lightDirection);

    Vec3f fuzzWeight = EvaluateSheen(
        viewDirection, lightDirection, bsdf.fuzzColor, bsdf.fuzzWeight, bsdf.fuzzRoughness);
    float fuzzPdf = PDFSheen(viewDirection, lightDirection, bsdf.fuzzRoughness);

    float totalPdf = substratePdf * (1.f - sheenProb) + fuzzPdf * sheenProb;

    // NOTE: scale albedo of base substrate to account for the light reflected by the fuzz layer
    float substrateScale = 1.f - fuzzReflectance * bsdf.fuzzWeight;
    weight = (substrateScale * substrateWeight + fuzzWeight) / totalPdf;

    // if (sheenProb > 0.f)
    // {
    //     printf("%f %f %f, %f %f %f, %f %f, %f, %f\n",
    //            substrateWeight.x,
    //            substrateWeight.y,
    //            substrateWeight.z,
    //            fuzzWeight.x,
    //            fuzzWeight.y,
    //            fuzzWeight.z,
    //            substratePdf,
    //            fuzzPdf,
    //            substrateScale,
    //            totalPdf);
    // }

    return true;
}

YBI_GPU Vec3f EvaluateOpenPBR(const OpenPbrBsdf &bsdf,
                              const Vec3f &viewDirection,
                              const Vec3f &lightDirection)
{
    Vec3f ltcCoefficients = GetLTCCoefficients(bsdf.fuzzRoughness, viewDirection.z);
    float fuzzReflectance = ltcCoefficients[2];

    Vec3f substrateWeight = EvaluateOpenPBRBaseSubstrate(bsdf, viewDirection, lightDirection);
    Vec3f fuzzWeight = EvaluateSheen(
        viewDirection, lightDirection, bsdf.fuzzColor, bsdf.fuzzWeight, bsdf.fuzzRoughness);

    // NOTE: scale albedo of base substrate to account for the light reflected by the fuzz layer
    float substrateScale = 1.f - fuzzReflectance * bsdf.fuzzWeight;
    Vec3f weight = substrateScale * substrateWeight + fuzzWeight;
    return weight;
}

YBI_GPU float
PDFOpenPBR(const OpenPbrBsdf &bsdf, const Vec3f &viewDirection, const Vec3f &lightDirection)
{
    float substratePdf = PDFOpenPBRBaseSubstrate(bsdf, viewDirection, lightDirection);
    float fuzzPdf = PDFSheen(viewDirection, lightDirection, bsdf.fuzzRoughness);

    // TODO: stop copy pasting this
    Vec3f ltcCoefficients = GetLTCCoefficients(bsdf.fuzzRoughness, viewDirection.z);
    float fuzzReflectance = ltcCoefficients[2];

    float sheenWeight = bsdf.fuzzWeight * fuzzReflectance * bsdf.fuzzColor.MaxComponent();
    float sheenProb = sheenWeight / (1.f + sheenWeight);

    float totalPdf = substratePdf * (1.f - sheenProb) + fuzzPdf * sheenProb;
    return totalPdf;
}

YBI_GPU bool SampleOpenPBRBaseSubstrate(const OpenPbrBsdf &bsdf,
                                        const Vec3f &viewDirection,
                                        const Vec3f &throughput,
                                        Vec2f u,
                                        Vec3f &lightDirection,
                                        Vec3f &weight)
{
    bool success = true;
    int lobeIndex = 0;
    if (u.x < 0.5f)
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
}

YBI_GPU Vec3f EvaluateOpenPBRBaseSubstrate(const OpenPbrBsdf &bsdf,
                                           const Vec3f &viewDirection,
                                           const Vec3f &lightDirection)
{
    Vec3f weight(0.f);
    // TODO: not hardcoded 0. also the diffuse lobe should be attenuated by gloss
    // reflection
    weight += EvaluateDiffuseLobe(viewDirection, lightDirection, bsdf.baseColor, 0.f);
    weight += EvaluateGlossLobe(viewDirection,
                                lightDirection,
                                bsdf.specularColor,
                                bsdf.specularIor,
                                bsdf.specularRoughness);
    return weight;
}

YBI_GPU float PDFOpenPBRBaseSubstrate(const OpenPbrBsdf &bsdf,
                                      const Vec3f &viewDirection,
                                      const Vec3f &lightDirection)
{
    float diffusePdf = OrenNayarPDF(lightDirection);
    float specularPdf = PDFGGX(viewDirection, lightDirection, bsdf.specularRoughness);

    float pdf = diffusePdf * 0.5f + specularPdf * 0.5f;

    pdf = Max(0.f, pdf);

    return pdf;
}

} // namespace ybi
