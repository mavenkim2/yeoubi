#pragma once

#include "core/frame.h"
#include "render/shaders/shading.h"

namespace ybi
{

YBI_GPU OrthonormalBasis MaterialXBsdfFrame(const ShadingInfo &info,
                                            Vec3f normal,
                                            const Vec3f &worldViewDirection)
{
    normal = LengthSquared(normal) > 1e-12f ? Normalize(normal) : info.shadingNormal;

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

    if (Dot(worldViewDirection, normal) < 0.f)
    {
        normal = -normal;
        bitangent = -bitangent;
    }
    return OrthonormalBasis(tangent, bitangent, normal);
}

} // namespace ybi
