#pragma once

#include "render/bsdf/openpbr/bsdf_types.h"
#include "render/bsdf/openpbr/ggx.h"
#include "render/bsdf/openpbr/sample_table.h"
#include "render/bsdf/openpbr/substrate.h"

namespace ybi
{

YBI_GPU float CalculateCoatEnergyReflected(
    uint64_t texture, float ior, float cosTheta, float isotropicAlpha, float coatWeight)
{
    // Get precomputed directional reflectance from table
    float energyReflected = 0.f;
    if (cosTheta > 0.f)
    {
        if (isotropicAlpha == 0.f)
        {
            energyReflected = FrDielectric(cosTheta, ior);
        }
        else
        {
            energyReflected =
                1.f - SampleOpaqueDielectricTable(texture, ior, cosTheta, Sqr(isotropicAlpha));
        }
    }

    energyReflected *= coatWeight;
    return energyReflected;
}

YBI_GPU float CalculateCoatLobeProbability(const OpenPbrBsdf &bsdf)
{
    float baseLobeContribution = 0.f;
    for (int i = 0; i < LobeType::COUNT; i++)
    {
        baseLobeContribution += bsdf.lobeWeights[i];
    }

    baseLobeContribution *= bsdf.coatDarkeningFactor.MaxComponent();

    float sum = bsdf.coatEnergyReflected + baseLobeContribution;
    float pdf = sum == 0.f ? 0.5 : bsdf.coatEnergyReflected / sum;
    return pdf;
}

// NOTE: See 3.4.2 View-dependent absorption
YBI_GPU Vec3f ScaleCoatColor(const Vec3f &coatColor, float cosTheta, float ior, float coatWeight)
{
    if (cosTheta > 0.f && coatColor.MinComponent() < 1.0f)
    {
        const Vec3f transmissionAtNormalIncidence = Sqrt(coatColor);
        const float refractedCosine = sqrtf(Max(0.f, 1.f - (1.f - Sqr(cosTheta)) / Sqr(ior)));
        const float distanceScale = 1.f / refractedCosine;
        const Vec3f coatTransmissionAlongRefractedRay =
            Pow(transmissionAtNormalIncidence, distanceScale);
        return Lerp(1.f, coatTransmissionAlongRefractedRay, coatWeight);
    }
    else
    {
        return Vec3f(1.f);
    }
}

YBI_GPU bool SampleCoatLayer(const OpenPbrBsdf &bsdf,
                             const Vec3f &viewDirection,
                             const Vec3f &throughput,
                             Vec2f u,
                             Vec3f &lightDirection,
                             Vec3f &weight)
{
    float coatProb = CalculateCoatLobeProbability(bsdf);

    if (bsdf.inside)
    {
        bool success =
            SampleOpenPBRBaseSubstrate(bsdf, viewDirection, throughput, u, lightDirection, weight);

        if (success)
        {
            const Vec3f baseLayerScale =
                ScaleCoatColor(bsdf.coatColor, -lightDirection.z, bsdf.coatIor, bsdf.coatWeight);
            weight *= baseLayerScale;
        }
        return success;
    }
    else
    {
        bool success = true;
        if (u.x < coatProb)
        {
            u.x = u.x / coatProb;
            u.x = ClampRandomNumber(u.x);

            success = SampleGGX(
                viewDirection, bsdf.coatRoughness, u, bsdf.coatIor, lightDirection, weight);
        }
        else
        {
            u.x = (u.x - coatProb) / (1 - coatProb);
            u.x = ClampRandomNumber(u.x);
            success = SampleOpenPBRBaseSubstrate(
                bsdf, viewDirection, throughput, u, lightDirection, weight);
        }

        if (!success)
            return false;

        // Evaluate coat and base layers
        Vec3f coatWeight =
            EvaluateGGX(viewDirection, lightDirection, bsdf.coatIor, bsdf.coatRoughness);
        float coatPdf = PDFGGX(viewDirection, lightDirection, bsdf.coatRoughness);

        Vec3f baseWeight = EvaluateOpenPBRBaseSubstrate(bsdf, viewDirection, lightDirection);
        float basePdf = PDFOpenPBRBaseSubstrate(bsdf, viewDirection, lightDirection);

        Vec3f totalWeight =
            coatWeight * bsdf.coatWeight +
            baseWeight * bsdf.coatDarkeningFactor *
                ScaleCoatColor(bsdf.coatColor, lightDirection.z, bsdf.coatIor, bsdf.coatWeight);

        float totalPdf = Lerp(basePdf, coatPdf, coatProb);

        weight = totalWeight / totalPdf;
        return true;
    }
}

YBI_GPU Vec3f EvaluateCoatLayer(const OpenPbrBsdf &bsdf,
                                const Vec3f &viewDirection,
                                const Vec3f &lightDirection)
{
    Vec3f baseWeight = EvaluateOpenPBRBaseSubstrate(bsdf, viewDirection, lightDirection);

    Vec3f totalWeight;
    if (bsdf.inside)
    {
        totalWeight =
            baseWeight *
            ScaleCoatColor(bsdf.coatColor, -lightDirection.z, bsdf.coatIor, bsdf.coatWeight);
    }
    else
    {
        Vec3f coatWeight =
            EvaluateGGX(viewDirection, lightDirection, bsdf.coatIor, bsdf.coatRoughness);

        totalWeight =
            coatWeight * bsdf.coatWeight +
            baseWeight * bsdf.coatDarkeningFactor *
                ScaleCoatColor(bsdf.coatColor, lightDirection.z, bsdf.coatIor, bsdf.coatWeight);
    }
    return totalWeight;
}

YBI_GPU float
PDFCoatLayer(const OpenPbrBsdf &bsdf, const Vec3f &viewDirection, const Vec3f &lightDirection)
{
    float totalPdf;
    float basePdf = PDFOpenPBRBaseSubstrate(bsdf, viewDirection, lightDirection);

    if (bsdf.inside)
    {
        totalPdf = basePdf;
    }
    else
    {
        float coatPdf = PDFGGX(viewDirection, lightDirection, bsdf.coatRoughness);
        float coatProb = CalculateCoatLobeProbability(bsdf);
        totalPdf = Lerp(basePdf, coatPdf, coatProb);
    }
    return totalPdf;
}

} // namespace ybi
