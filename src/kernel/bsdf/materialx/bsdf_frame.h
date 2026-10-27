#pragma once

#include "core/frame.h"
#include "kernel/bsdf/closure_types.h"

namespace ybi
{

YBI_GPU OrthonormalBasis MaterialXBsdfFrame(const ClosureEvaluationInfo &info,
                                            Vec3f normal,
                                            const Vec3f &worldViewDirection,
                                            bool *frontFacing = nullptr)
{
    normal = LengthSquared(normal) > 1e-12f ? Normalize(normal) : info.shadingNormal;
    if (Dot(normal, info.geometricNormal) < 0.f)
        normal = -normal;

    Vec3f tangent;
    Vec3f bitangent = Cross(normal, info.dpdu);
    if (LengthSquared(bitangent) > 1e-12f)
    {
        bitangent = Normalize(bitangent);
        tangent = Normalize(Cross(bitangent, normal));
    }
    else
    {
        BuildOrthonormalBasis(normal, tangent, bitangent);
    }

    const float nDotV = Dot(worldViewDirection, normal);
    if (frontFacing)
        *frontFacing = nDotV > 0.f;
    if (nDotV < 0.f)
    {
        normal = -normal;
        bitangent = -bitangent;
    }

    return OrthonormalBasis(tangent, bitangent, normal);
}

} // namespace ybi
