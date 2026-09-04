#include "core/assert.h"
#include "core/math.h"
#include "core/vec2f.h"
#include "core/vec3f.h"
#include "device/cuda_compat.h"
#include "render/bsdf/ltc_coefficients.h"
#include "render/sampling.h"
#include <cmath>

namespace ybi
{

YBI_GPU float ClampRandomNumber(float u)
{
    YBI_ASSERT(u >= 0.f);
    u = Min(u, FLOAT_BEFORE_ONE);
    return u;
}

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

YBI_CPU_GPU Vec3f SampleSheenLTC(Vec2f u, float roughness)
{
    Vec3f w = SampleCosineHemisphere(u);
    Vec3f ltcCoeffs = GetLTCCoefficients(roughness, w.z);

    float aInv = ltcCoeffs[0];
    float bInv = ltcCoeffs[1];
    float a = 1.f / aInv;

    Vec3f result(a * w.x - bInv * a * w.z, a * w.y, w.z);
    return result;
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

YBI_CPU_GPU Vec3f EvaluateSheen(const Vec3f &wo, const Vec3f &wi, float roughness)
{
    if (wo.z <= 0.f || wi.z <= 0.f)
    {
        return Vec3f(0.f);
    }

    Vec3f wiStd = RotateVector(wi, Vec2f(wo.x, -wo.y));
    Vec3f ltcCoeffs = GetLTCCoefficients(roughness, wo.z);

    Vec3f value(EvaluateLTC(wiStd, ltcCoeffs));

    // TODO: get lobe
    // value *= fuzzWeight * fuzzColor * ltcCoeffs[2];
    return value;
}

YBI_CPU_GPU float PDFSheen(const Vec3f &wo, const Vec3f &wi)
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

YBI_CPU_GPU void SampleOpenPBRFuzz(const Vec3f &wo, Vec2f u)
{
    // TODO: calculate the probability of sampling from this fuzz layer
    float prob = 0.f;

    if (u.x < prob)
    {
        u.x /= prob;
        u.x = ClampRandomNumber(u.x);

        Vec3f wi = SampleSheenLTC(u);
        wi = RotateVector(wi, wo.Xy());
    }
    else
    {
    }
}

YBI_CPU_GPU Vec3f SampleGGXVNDF(const Vec3f &wo, const Vec2f &a, const Vec2f &u)
{
    // fma(1 - u.y, 1 + wo.z, -wo.z)
    // z = 1 + wo.z - u.y - u.y * wo.z - wo.z = 1 - u.y - u.y * wo.z

    Vec3f woStd = Normalize(Vec3f(a * wo.Xy(), wo.z));

    float z = 1 - u.y - u.y * woStd.z;
    float sinTheta = Clamp(sqrtf(1 - z * z), 0.f, 1.f);

    float phi = 2 * YBI_PI * u.x;
    float cosPhi = Cos(phi);
    float sinPhi = Sin(phi);

    Vec3f result(sinTheta * cosPhi, sinTheta * sinPhi, z);
    result += woStd;

    result.x *= a.x;
    result.y *= a.y;
    return Normalize(result);
}

YBI_CPU_GPU float EvaluateMaskingGGX() {}

YBI_CPU_GPU float EvaluateAnisoGGX(const Vec3f &v, const Vec3f &a)
{
    float x = (v.x / a.x) * (v.x / a.x);
    float y = (v.y / a.y) * (v.y / a.y);
    float z = v.z * v.z;
    float denom = x + y + z;
    denom *= denom;

    float result = 1.f / (YBI_PI * a.x * a.y * denom);
    return result;
}

YBI_CPU_GPU float FrDielectric(float eta) {}

YBI_CPU_GPU bool
SampleGGX(const Vec3f &viewDirection, const Vec2f &a, const Vec2f &u, Vec3f &lightDirection)
{
    Vec3f halfway = SampleGGXVNDF(viewDirection, a, u);
    float hDotV = Dot(halfway, viewDirection);

    if (hDotV < 0.f)
    {
        return false;
    }

    lightDirection = Normalize(-viewDirection + 2.f * hDotV * halfway);
    float lDotV = Dot(lightDirection, viewDirection);

    if (lDotV < 0.f)
    {
        return false;
    }

    float D = EvaluateAnisoGGX(wm, a);
    // float G = ;
}

} // namespace ybi
