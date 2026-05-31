#pragma once

#include "render/integrator_common.h"
#include "render/launch_params.h"

namespace ybi
{
namespace render
{
namespace integrator
{

struct RayDifferential
{
    Vec3 origin = Vec3(0.0f, 0.0f, 0.0f);
    Vec3 dir = Vec3(0.0f, 0.0f, 1.0f);
    Vec3 originX = Vec3(0.0f, 0.0f, 0.0f);
    Vec3 originY = Vec3(0.0f, 0.0f, 0.0f);
    Vec3 dirX = Vec3(0.0f, 0.0f, 1.0f);
    Vec3 dirY = Vec3(0.0f, 0.0f, 1.0f);
    bool valid = false;
};

struct TextureDifferentialResult
{
    float dSdx = 0.0f;
    float dTdx = 0.0f;
    float dSdy = 0.0f;
    float dTdy = 0.0f;
    bool valid = false;
};

struct HitPlaneDifferentialResult
{
    Vec3 dpdx = Vec3(0.0f, 0.0f, 0.0f);
    Vec3 dpdy = Vec3(0.0f, 0.0f, 0.0f);
    bool valid = false;
};

YBI_DEVICE bool IsFiniteScalar(float value)
{
    return value == value && fabsf(value) <= 1.0e30f;
}

YBI_DEVICE bool IsFiniteVec3(const Vec3 &value)
{
    return IsFiniteScalar(value.x) && IsFiniteScalar(value.y) && IsFiniteScalar(value.z);
}

YBI_DEVICE Float4x4 RotateFromTo(const Vec3 &from, const Vec3 &to)
{
    if (!IsFiniteVec3(from) || !IsFiniteVec3(to) || LengthSquared(from) <= 1.0e-20f ||
        LengthSquared(to) <= 1.0e-20f)
    {
        return Float4x4::Identity();
    }

    const Vec3 fromN = Normalize(from);
    const Vec3 toN = Normalize(to);
    if (Dot(fromN, toN) >= 1.0f - 1.0e-6f)
    {
        return Float4x4::Identity();
    }

    Vec3 refl = Vec3(0.0f, 0.0f, 1.0f);
    if (fabsf(fromN.x) < 0.72f && fabsf(toN.x) < 0.72f)
    {
        refl = Vec3(1.0f, 0.0f, 0.0f);
    }
    else if (fabsf(fromN.y) < 0.72f && fabsf(toN.y) < 0.72f)
    {
        refl = Vec3(0.0f, 1.0f, 0.0f);
    }

    const Vec3 u = refl - fromN;
    const Vec3 v = refl - toN;
    const float uu = Dot(u, u);
    const float vv = Dot(v, v);
    if (uu <= 1.0e-20f || vv <= 1.0e-20f)
    {
        return Float4x4::Identity();
    }

    const float uv = Dot(u, v);
    const float uuScale = 2.0f / uu;
    const float vvScale = 2.0f / vv;
    const float uvScale = 4.0f * uv / (uu * vv);
    Float4x4 r = Float4x4::Identity();
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            r.m[i][j] = (i == j ? 1.0f : 0.0f) - uuScale * u[i] * u[j] - vvScale * v[i] * v[j] +
                        uvScale * v[i] * u[j];
        }
    }
    return r;
}

YBI_DEVICE void InvalidateRayDifferential(RayDifferential *rayDiff)
{
    if (!rayDiff)
    {
        return;
    }

    rayDiff->originX = rayDiff->origin;
    rayDiff->originY = rayDiff->origin;
    rayDiff->dirX = rayDiff->dir;
    rayDiff->dirY = rayDiff->dir;
    rayDiff->valid = false;
}

YBI_DEVICE bool SolveLinear2x2(float a00,
                                      float a01,
                                      float a10,
                                      float a11,
                                      float b0,
                                      float b1,
                                      float *outX,
                                      float *outY)
{
    if (!outX || !outY)
    {
        return false;
    }

    const float det = a00 * a11 - a01 * a10;
    if (fabsf(det) <= 1.0e-8f)
    {
        return false;
    }

    const float invDet = 1.0f / det;
    const float x = (b0 * a11 - a01 * b1) * invDet;
    const float y = (a00 * b1 - b0 * a10) * invDet;
    if (!IsFiniteScalar(x) || !IsFiniteScalar(y))
    {
        return false;
    }

    *outX = x;
    *outY = y;
    return true;
}

YBI_DEVICE TextureDifferentialResult ComputeTextureDifferentials(const Vec3 &dpdx,
                                                                        const Vec3 &dpdy,
                                                                        const Vec3 &dPds,
                                                                        const Vec3 &dPdt,
                                                                        const Vec3 &geomNormal)
{
    TextureDifferentialResult result = {};
    if (!IsFiniteVec3(dpdx) || !IsFiniteVec3(dpdy) || !IsFiniteVec3(dPds) ||
        !IsFiniteVec3(dPdt) || !IsFiniteVec3(geomNormal))
    {
        return result;
    }

    const float dPdsLengthSq = Dot(dPds, dPds);
    const float dPdtLengthSq = Dot(dPdt, dPdt);
    const float geomNormalLengthSq = Dot(geomNormal, geomNormal);
    if (dPdsLengthSq <= 1.0e-8f || dPdtLengthSq <= 1.0e-8f || geomNormalLengthSq <= 1.0e-8f)
    {
        return result;
    }

    const float basisDot = Dot(dPds, dPdt);
    if (fabsf(basisDot) <= 1.0e-8f)
    {
        const float invDsds = 1.0f / dPdsLengthSq;
        const float invDtdt = 1.0f / dPdtLengthSq;
        result.dSdx = Dot(dpdx, dPds) * invDsds;
        result.dTdx = Dot(dpdx, dPdt) * invDtdt;
        result.dSdy = Dot(dpdy, dPds) * invDsds;
        result.dTdy = Dot(dpdy, dPdt) * invDtdt;
    }
    else
    {
        const Vec3 n2 = geomNormal * geomNormal;
        float a00 = 0.0f;
        float a01 = 0.0f;
        float a10 = 0.0f;
        float a11 = 0.0f;
        float bx0 = 0.0f;
        float bx1 = 0.0f;
        float by0 = 0.0f;
        float by1 = 0.0f;

        if (n2.x > n2.y && n2.x > n2.z)
        {
            a00 = dPds.y;
            a01 = dPdt.y;
            a10 = dPds.z;
            a11 = dPdt.z;
            bx0 = dpdx.y;
            bx1 = dpdx.z;
            by0 = dpdy.y;
            by1 = dpdy.z;
        }
        else if (n2.y > n2.z)
        {
            a00 = dPds.x;
            a01 = dPdt.x;
            a10 = dPds.z;
            a11 = dPdt.z;
            bx0 = dpdx.x;
            bx1 = dpdx.z;
            by0 = dpdy.x;
            by1 = dpdy.z;
        }
        else
        {
            a00 = dPds.x;
            a01 = dPdt.x;
            a10 = dPds.y;
            a11 = dPdt.y;
            bx0 = dpdx.x;
            bx1 = dpdx.y;
            by0 = dpdy.x;
            by1 = dpdy.y;
        }

        if (!SolveLinear2x2(
                a00, a01, a10, a11, bx0, bx1, &result.dSdx, &result.dTdx) ||
            !SolveLinear2x2(
                a00, a01, a10, a11, by0, by1, &result.dSdy, &result.dTdy))
        {
            result.dSdx = 0.0f;
            result.dTdx = 0.0f;
            result.dSdy = 0.0f;
            result.dTdy = 0.0f;
            return result;
        }
    }

    result.valid = IsFiniteScalar(result.dSdx) && IsFiniteScalar(result.dTdx) &&
                   IsFiniteScalar(result.dSdy) && IsFiniteScalar(result.dTdy);
    if (!result.valid)
    {
        result.dSdx = 0.0f;
        result.dTdx = 0.0f;
        result.dSdy = 0.0f;
        result.dTdy = 0.0f;
    }
    return result;
}

YBI_DEVICE bool ApproximateDpDxy(const LaunchParams &params,
                                        const Vec3 &p,
                                        const Vec3 &n,
                                        int samplesPerPixel,
                                        Vec3 *outDpdx,
                                        Vec3 *outDpdy)
{
    if (!outDpdx || !outDpdy || !IsFiniteVec3(p) || !IsFiniteVec3(n) ||
        !IsFiniteVec3(params.minPosDifferentialX) || !IsFiniteVec3(params.minPosDifferentialY) ||
        !IsFiniteVec3(params.minDirDifferentialX) || !IsFiniteVec3(params.minDirDifferentialY))
    {
        return false;
    }

    *outDpdx = Vec3(0.0f);
    *outDpdy = Vec3(0.0f);

    const Vec3 pCamera = TransformPointAffine(params.cameraFromWorld, p);
    const Vec3 nCamera = TransformVectorAffine(params.cameraFromWorld, n);
    if (!IsFiniteVec3(pCamera) || !IsFiniteVec3(nCamera) || LengthSquared(pCamera) <= 1.0e-20f ||
        LengthSquared(nCamera) <= 1.0e-20f)
    {
        return false;
    }

    const Float4x4 downZFromCamera = RotateFromTo(Normalize(pCamera), Vec3(0.0f, 0.0f, 1.0f));
    const Vec3 pDownZ = TransformPointAffine(downZFromCamera, pCamera);
    const Vec3 nDownZ = TransformVectorAffine(downZFromCamera, nCamera);
    const float d = Dot(nDownZ, pDownZ);

    const Vec3 xRayOrigin = params.minPosDifferentialX;
    const Vec3 yRayOrigin = params.minPosDifferentialY;
    const Vec3 xRayDir = Vec3(0.0f, 0.0f, 1.0f) + params.minDirDifferentialX;
    const Vec3 yRayDir = Vec3(0.0f, 0.0f, 1.0f) + params.minDirDifferentialY;
    const float xDenom = Dot(nDownZ, xRayDir);
    const float yDenom = Dot(nDownZ, yRayDir);
    if (fabsf(xDenom) <= 1.0e-8f || fabsf(yDenom) <= 1.0e-8f)
    {
        return false;
    }

    const float tx = -(Dot(nDownZ, xRayOrigin) - d) / xDenom;
    const float ty = -(Dot(nDownZ, yRayOrigin) - d) / yDenom;
    if (!IsFiniteScalar(tx) || !IsFiniteScalar(ty))
    {
        return false;
    }

    const Vec3 px = xRayOrigin + tx * xRayDir;
    const Vec3 py = yRayOrigin + ty * yRayDir;
    const Float4x4 cameraFromDownZ = Transpose(downZFromCamera);
    const Vec3 dpdxCamera = TransformVectorAffine(cameraFromDownZ, px - pDownZ);
    const Vec3 dpdyCamera = TransformVectorAffine(cameraFromDownZ, py - pDownZ);
    const float sppScale = fmaxf(0.125f, 1.0f / sqrtf(float(Max(samplesPerPixel, 1))));
    *outDpdx = sppScale * TransformVectorAffine(params.worldFromCamera, dpdxCamera);
    *outDpdy = sppScale * TransformVectorAffine(params.worldFromCamera, dpdyCamera);
    return IsFiniteVec3(*outDpdx) && IsFiniteVec3(*outDpdy);
}

YBI_DEVICE bool TransferRayDifferentialToHitPoint(const RayDifferential &rayDiff,
                                                         const Vec3 &hitPoint,
                                                         const Vec3 &geomNormal,
                                                         Vec3 *outPx,
                                                         Vec3 *outPy)
{
    if (!rayDiff.valid || !outPx || !outPy || !IsFiniteVec3(hitPoint) || !IsFiniteVec3(geomNormal) ||
        !IsFiniteVec3(rayDiff.originX) || !IsFiniteVec3(rayDiff.originY) ||
        !IsFiniteVec3(rayDiff.dirX) || !IsFiniteVec3(rayDiff.dirY) ||
        LengthSquared(geomNormal) <= 1.0e-20f)
    {
        return false;
    }

    const float dxDenom = Dot(geomNormal, rayDiff.dirX);
    const float dyDenom = Dot(geomNormal, rayDiff.dirY);
    if (fabsf(dxDenom) <= 1.0e-8f || fabsf(dyDenom) <= 1.0e-8f)
    {
        return false;
    }

    const float tx = Dot(geomNormal, hitPoint - rayDiff.originX) / dxDenom;
    const float ty = Dot(geomNormal, hitPoint - rayDiff.originY) / dyDenom;
    if (!IsFiniteScalar(tx) || !IsFiniteScalar(ty))
    {
        return false;
    }

    *outPx = rayDiff.originX + rayDiff.dirX * tx;
    *outPy = rayDiff.originY + rayDiff.dirY * ty;
    return IsFiniteVec3(*outPx) && IsFiniteVec3(*outPy);
}

YBI_DEVICE HitPlaneDifferentialResult ComputeHitPlaneDifferentials(const LaunchParams &params,
                                                                          const RayDifferential &rayDiff,
                                                                          const Vec3 &hitPoint,
                                                                          const Vec3 &geomNormal,
                                                                          int samplesPerPixel)
{
    HitPlaneDifferentialResult result = {};
    if (!IsFiniteVec3(hitPoint) || !IsFiniteVec3(geomNormal))
    {
        return result;
    }

    Vec3 px = Vec3(0.0f);
    Vec3 py = Vec3(0.0f);
    if (TransferRayDifferentialToHitPoint(rayDiff, hitPoint, geomNormal, &px, &py))
    {
        result.dpdx = px - hitPoint;
        result.dpdy = py - hitPoint;
        result.valid = IsFiniteVec3(result.dpdx) && IsFiniteVec3(result.dpdy);
        return result;
    }

    if (ApproximateDpDxy(params, hitPoint, geomNormal, samplesPerPixel, &result.dpdx, &result.dpdy))
    {
        result.valid = IsFiniteVec3(result.dpdx) && IsFiniteVec3(result.dpdy);
    }
    return result;
}

YBI_DEVICE Vec3 ComputePerspectiveCameraDirection(const LaunchParams &params,
                                                         float rasterX,
                                                         float rasterY)
{
    const Vec3 cameraPoint =
        TransformPointPerspective(params.cameraFromRaster, Vec3(rasterX, rasterY, 0.0f));
    if (!(cameraPoint.x == cameraPoint.x && cameraPoint.y == cameraPoint.y &&
          cameraPoint.z == cameraPoint.z))
    {
        return Vec3(0.0f, 0.0f, 1.0f);
    }

    const Vec3 worldPoint = TransformPointAffine(params.worldFromCamera, cameraPoint);
    const Vec3 worldOrigin =
        Vec3(params.cameraOrigin.x, params.cameraOrigin.y, params.cameraOrigin.z);
    return Normalize(worldPoint - worldOrigin);
}

YBI_DEVICE RayDifferential InitPerspectiveRayDifferential(const LaunchParams &params,
                                                                 float pixelX,
                                                                 float pixelY,
                                                                 unsigned int width,
                                                                 unsigned int height)
{
    RayDifferential rayDiff = {};
    const Vec3 cameraOrigin(params.cameraOrigin.x, params.cameraOrigin.y, params.cameraOrigin.z);
    rayDiff.origin = cameraOrigin;
    rayDiff.originX = cameraOrigin;
    rayDiff.originY = cameraOrigin;

    if (width == 0u || height == 0u)
    {
        return rayDiff;
    }

    rayDiff.dir = ComputePerspectiveCameraDirection(params, pixelX + 0.5f, pixelY + 0.5f);
    rayDiff.dirX = ComputePerspectiveCameraDirection(params, pixelX + 1.5f, pixelY + 0.5f);
    rayDiff.dirY = ComputePerspectiveCameraDirection(params, pixelX + 0.5f, pixelY + 1.5f);
    rayDiff.valid = true;
    return rayDiff;
}

} // namespace integrator
} // namespace render
} // namespace ybi
