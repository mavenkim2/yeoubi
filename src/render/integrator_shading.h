#pragma once

#include "render/integrator_texture.h"

namespace ybi
{
namespace render
{
namespace integrator
{

struct ShadingFrame
{
    Vec3 geomNormal = Vec3(0.0f, 0.0f, 1.0f);
    Vec3 shadingNormal = Vec3(0.0f, 0.0f, 1.0f);
    Vec3 tangent = Vec3(1.0f, 0.0f, 0.0f);
    Vec3 bitangent = Vec3(0.0f, 1.0f, 0.0f);
    float handedness = 1.0f;
    bool hasTangentFrame = false;
};

YBI_INTEGRATOR_HD Vec3 SafeNormalizeOrDefault(const Vec3 &value, const Vec3 &fallback)
{
    return LengthSquared(value) > 1.0e-12f ? Normalize(value) : fallback;
}

YBI_INTEGRATOR_HD bool TryComputeTriangleTangentFrame(
    const LaunchParams::InstanceGeomRef &geomRef,
    const HitInfo &hit,
    const Vec3 &baseNormal,
    Vec3 *outTangent,
    Vec3 *outBitangent,
    float *outHandedness)
{
    if (!outTangent || !outBitangent || !outHandedness || !hit.hasWorldTriangle || !hit.hasBarycentrics)
    {
        return false;
    }

    MaterialTextureSampleInputs inputs = {};
    if (!TryComputeMaterialTextureSampleInputs(
            geomRef, static_cast<unsigned int>(hit.primitiveIndex), hit.barycentrics, &inputs))
    {
        return false;
    }

    const Vec3 edge1 = hit.worldTri1 - hit.worldTri0;
    const Vec3 edge2 = hit.worldTri2 - hit.worldTri0;
    const Vec2 duv1(inputs.uv1.x - inputs.uv0.x, inputs.uv1.y - inputs.uv0.y);
    const Vec2 duv2(inputs.uv2.x - inputs.uv0.x, inputs.uv2.y - inputs.uv0.y);
    const float det = duv1.x * duv2.y - duv1.y * duv2.x;
    if (fabsf(det) <= 1.0e-8f)
    {
        return false;
    }

    const float invDet = 1.0f / det;
    Vec3 tangent = (edge1 * duv2.y - edge2 * duv1.y) * invDet;
    Vec3 bitangent = (edge2 * duv1.x - edge1 * duv2.x) * invDet;
    tangent = tangent - baseNormal * Dot(tangent, baseNormal);
    bitangent = bitangent - baseNormal * Dot(bitangent, baseNormal);
    if (LengthSquared(tangent) <= 1.0e-12f || LengthSquared(bitangent) <= 1.0e-12f)
    {
        return false;
    }

    tangent = Normalize(tangent);
    bitangent = Normalize(bitangent);
    *outHandedness = Dot(Cross(tangent, bitangent), baseNormal) >= 0.0f ? 1.0f : -1.0f;
    *outTangent = tangent;
    *outBitangent = Normalize(Cross(baseNormal, tangent)) * *outHandedness;
    return true;
}

template <typename State>
YBI_INTEGRATOR_HD ShadingFrame EvaluateShadingFrame(State &state,
                                                    const LaunchParams::InstanceGeomRef &geomRef,
                                                    const HitInfo &hit)
{
    ShadingFrame frame = {};
    frame.geomNormal = hit.hasGeomNormal
                           ? SafeNormalizeOrDefault(hit.geomNormal, Vec3(0.0f, 0.0f, 1.0f))
                           : SafeNormalizeOrDefault(-hit.rayDir, Vec3(0.0f, 0.0f, 1.0f));

    Vec3 baseNormal =
        hit.hasShadingNormal ? SafeNormalizeOrDefault(hit.shadingNormal, frame.geomNormal) : frame.geomNormal;
    if (Dot(baseNormal, frame.geomNormal) < 0.0f)
    {
        baseNormal = -baseNormal;
    }

    if (!TryComputeTriangleTangentFrame(
            geomRef, hit, baseNormal, &frame.tangent, &frame.bitangent, &frame.handedness))
    {
        BuildOrthonormalBasis(baseNormal, frame.tangent, frame.bitangent);
        frame.handedness = 1.0f;
        frame.hasTangentFrame = false;
    }
    else
    {
        frame.hasTangentFrame = true;
    }

    frame.shadingNormal = baseNormal;

    Vec4 normalSample = {};
    if (frame.hasTangentFrame &&
        TrySampleMaterialTextureSemantic(state, geomRef, hit, kSemanticNormal, normalSample))
    {
        if (normalSample.x == 0.0f && normalSample.y == 0.0f && normalSample.z == 0.0f &&
            normalSample.w == 0.0f)
        {
            return frame;
        }
        Vec3 tangentNormal(
            normalSample.x * 2.0f - 1.0f, normalSample.y * 2.0f - 1.0f, normalSample.z * 2.0f - 1.0f);
        if (LengthSquared(tangentNormal) > 1.0e-12f)
        {
            tangentNormal = Normalize(tangentNormal);
            const Vec3 mappedNormal = Normalize(frame.tangent * tangentNormal.x +
                                                frame.bitangent * tangentNormal.y +
                                                baseNormal * tangentNormal.z);
            if (Dot(mappedNormal, frame.geomNormal) > 0.0f)
            {
                frame.shadingNormal = mappedNormal;
                frame.tangent = Normalize(frame.tangent -
                                          frame.shadingNormal * Dot(frame.tangent, frame.shadingNormal));
                frame.bitangent = Normalize(Cross(frame.shadingNormal, frame.tangent)) * frame.handedness;
            }
        }
    }

    return frame;
}

} // namespace integrator
} // namespace render
} // namespace ybi
