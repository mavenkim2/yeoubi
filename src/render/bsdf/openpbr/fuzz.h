#pragma once

#include "core/assert.h"
#include "core/math.h"
#include "core/vec2f.h"
#include "core/vec3f.h"
#include "device/cuda_compat.h"
#include "render/bsdf/ltc_coefficients.h"
#include "render/bsdf/openpbr/coat.h"
#include "render/sampling.h"

namespace ybi
{

YBI_GPU Vec3f RotateVector(const Vec3f &v, const Vec2f &woXY)
{
    const float lengthSquared = Dot(woXY, woXY);
    if (lengthSquared == 0.0f)
    {
        return v;
    }
    const float invLength = RcpSqrt(lengthSquared);
    const float sinPhi = woXY.y * invLength;
    const float cosPhi = woXY.x * invLength;

    return Vec3f(cosPhi * v.x + sinPhi * -v.y, cosPhi * v.y + sinPhi * v.x, v.z);
}

YBI_GPU Vec3f GetLTCCoefficients(float roughness, float cosTheta)
{
    // TODO: GPU hardware textures
    float row = Clamp(roughness, 0.f, FLOAT_BEFORE_ONE) * (LTC_ARRAY_DIM - 1);
    float col = Clamp(cosTheta, 0.f, FLOAT_BEFORE_ONE) * (LTC_ARRAY_DIM - 1);

    int rowIndex = static_cast<int>(row);
    int colIndex = static_cast<int>(col);

    float rowWeight = row - Floor(row);
    float colWeight = col - Floor(col);

    Vec3f r0 = Lerp(LTC_ARRAY[LTC_ARRAY_DIM * rowIndex + colIndex],
                    LTC_ARRAY[LTC_ARRAY_DIM * rowIndex + colIndex + 1],
                    colWeight);
    Vec3f r1 = Lerp(LTC_ARRAY[LTC_ARRAY_DIM * (rowIndex + 1) + colIndex],
                    LTC_ARRAY[LTC_ARRAY_DIM * (rowIndex + 1) + colIndex + 1],
                    colWeight);

    Vec3f result = Lerp(r0, r1, rowWeight);
    return result;
}

YBI_CPU_GPU Vec3f SampleSheenLTC(const Vec3f &ltcCoeffs, Vec2f u)
{
    Vec3f w = SampleCosineHemisphere(u);

    float aInv = ltcCoeffs[0];
    float bInv = ltcCoeffs[1];
    float a = 1.f / aInv;

    Vec3f result(a * w.x - bInv * a * w.z, a * w.y, w.z);
    return Normalize(result);
}

YBI_CPU_GPU float EvaluateLTC(const Vec3f &wiStd, const Vec3f &ltcCoeffs)
{
    float aInv = ltcCoeffs[0];
    float bInv = ltcCoeffs[1];

    Vec3f wiLtc(aInv * wiStd.x + bInv * wiStd.z, aInv * wiStd.y, wiStd.z);

    float z = Max(0.f, wiLtc.z);
    float lenSquared = Dot(wiLtc, wiLtc);
    if (lenSquared == 0.f || aInv == 0.f || z == 0.f)
        return 0.f;

    float aInvOverLenSquared = aInv / lenSquared;

    // Cos(wi) / PI / length^3
    return YBI_RCP_PI * z * aInvOverLenSquared * aInvOverLenSquared;
}

YBI_CPU_GPU bool
SampleSheen(const Vec3f &viewDirection, const Vec2f &u, float roughness, Vec3f &lightDirection)
{
    if (viewDirection.z <= 0.f)
    {
        return false;
    }

    Vec3f ltcCoeffs = GetLTCCoefficients(roughness, viewDirection.z);
    Vec3f lightDirectionStd = SampleSheenLTC(ltcCoeffs, u);
    lightDirection = RotateVector(lightDirectionStd, viewDirection.Xy());

    if (viewDirection.z * lightDirection.z <= 0.f)
        return false;
    return true;
}

YBI_CPU_GPU Vec3f EvaluateSheen(
    const Vec3f &wo, const Vec3f &wi, const Vec3f &fuzzColor, float fuzzWeight, float roughness)
{
    if (wo.z <= 0.f || wi.z <= 0.f)
    {
        return Vec3f(0.f);
    }

    Vec3f wiStd = RotateVector(wi, Vec2f(wo.x, -wo.y));
    Vec3f ltcCoeffs = GetLTCCoefficients(roughness, wo.z);

    Vec3f value(EvaluateLTC(wiStd, ltcCoeffs));

    value *= fuzzWeight * fuzzColor * ltcCoeffs[2];
    return value;
}

YBI_CPU_GPU float PDFSheen(const Vec3f &wo, const Vec3f &wi, float roughness)
{
    if (wo.z <= 0.f || wi.z <= 0.f)
    {
        return 0.f;
    }

    Vec3f wiStd = RotateVector(wi, Vec2f(wo.x, -wo.y));
    Vec3f ltcCoeffs = GetLTCCoefficients(roughness, wo.z);

    float pdf = EvaluateLTC(wiStd, ltcCoeffs);
    return pdf;
}

YBI_GPU bool SampleFuzzLayer(const OpenPbrBsdf &bsdf,
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
    // Calculate specular lobe weight float Fr =
    // FrDielectric(viewDirection.z, bsdf.specularReflectionIor);
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

        bool success = SampleCoatLayer(bsdf, viewDirection, throughput, u, lightDirection, weight);
        // SampleOpenPBRBaseSubstrate(bsdf, viewDirection, throughput, u, lightDirection, weight);

        if (!success)
            return false;
    }

    Vec3f coatWeight = EvaluateOpenPBRBaseSubstrate(bsdf, viewDirection, lightDirection);
    float coatPdf = PDFCoatLayer(bsdf, viewDirection, lightDirection);

    Vec3f fuzzWeight = EvaluateSheen(
        viewDirection, lightDirection, bsdf.fuzzColor, bsdf.fuzzWeight, bsdf.fuzzRoughness);
    float fuzzPdf = PDFSheen(viewDirection, lightDirection, bsdf.fuzzRoughness);

    float totalPdf = coatPdf * (1.f - sheenProb) + fuzzPdf * sheenProb;

    // NOTE: scale albedo of base substrate to account for the light reflected by the fuzz layer
    float coatScale = 1.f - fuzzReflectance * bsdf.fuzzWeight;
    weight = (coatScale * coatWeight + fuzzWeight) / totalPdf;

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

YBI_GPU Vec3f EvaluateOpenPBRLocal(const OpenPbrBsdf &bsdf,
                                   const Vec3f &viewDirection,
                                   const Vec3f &lightDirection)
{
    Vec3f ltcCoefficients = GetLTCCoefficients(bsdf.fuzzRoughness, viewDirection.z);
    float fuzzReflectance = ltcCoefficients[2];

    Vec3f coatWeight = EvaluateCoatLayer(bsdf, viewDirection, lightDirection);
    Vec3f fuzzWeight = EvaluateSheen(
        viewDirection, lightDirection, bsdf.fuzzColor, bsdf.fuzzWeight, bsdf.fuzzRoughness);

    // NOTE: scale albedo of base substrate to account for the light reflected by the fuzz layer
    float coatScale = 1.f - fuzzReflectance * bsdf.fuzzWeight;
    Vec3f weight = coatScale * coatWeight + fuzzWeight;
    return weight;
}

YBI_GPU float PDFOpenPBRLocal(const OpenPbrBsdf &bsdf,
                              const Vec3f &viewDirection,
                              const Vec3f &lightDirection)
{
    float coatPdf = PDFCoatLayer(bsdf, viewDirection, lightDirection);
    float fuzzPdf = PDFSheen(viewDirection, lightDirection, bsdf.fuzzRoughness);

    // TODO: stop copy pasting this
    Vec3f ltcCoefficients = GetLTCCoefficients(bsdf.fuzzRoughness, viewDirection.z);
    float fuzzReflectance = ltcCoefficients[2];

    float sheenWeight = bsdf.fuzzWeight * fuzzReflectance * bsdf.fuzzColor.MaxComponent();
    float sheenProb = sheenWeight / (1.f + sheenWeight);

    float totalPdf = coatPdf * (1.f - sheenProb) + fuzzPdf * sheenProb;
    return totalPdf;
}

} // namespace ybi
