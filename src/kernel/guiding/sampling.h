#pragma once

#include "core/frame.h"
#include "kernel/guiding/kd_tree_traverse.h"
#include "kernel/guiding/vmf.h"
#include <cfloat>

namespace ybi
{

// Prepare a private mixture so that shading does not change the trained field.
YBI_GPU bool PrepareGuidingDistribution(const PathGuidingKdTreeParams &params,
                                        const Vec3f &position,
                                        VMMDistribution &distribution)
{
    const VMMDistribution *stored = TraverseKdTree(params, position);
    if (!stored)
        return false;

    distribution = *stored;
    YBI_ASSERT(distribution.numComponents <= VMM_MAX_NUM_COMPONENTS);
    const Vec3f shift = distribution.pivotPosition - position;
    for (uint32_t k = 0; k < distribution.numComponents; ++k)
    {
        const float distance = distribution.distance[k];
        if (!(distance > 0.f) || !IsFinite(distance))
            continue;

        const Vec3f offset = distribution.meanDirection[k] * distance + shift;
        const float shiftedDistance = Length(offset);
        if (shiftedDistance > FLT_EPSILON && IsFinite(shiftedDistance))
        {
            distribution.meanDirection[k] = offset / shiftedDistance;
            distribution.distance[k] = shiftedDistance;
        }
    }
    distribution.pivotPosition = position;
    return true;
}

YBI_GPU Vec3f SampleVMF(const Vec3f &meanDirection, float kappa, const Vec2f &u)
{
    YBI_ASSERT(kappa >= 0.f && IsFinite(kappa));
    const float random = Clamp(u.x, 0.f, FLOAT_BEFORE_ONE);
    // expm1/log1p preserve precision for broad lobes with small kappa.
    const float cosTheta = Clamp(
        kappa == 0.f ? 1.f - 2.f * random : 1.f + log1pf(random * expm1f(-2.f * kappa)) / kappa,
        -1.f,
        1.f);
    const float sinTheta = Sqrt(Max(0.f, (1.f - cosTheta) * (1.f + cosTheta)));
    const float phi = 2.f * YBI_PI * Clamp(u.y, 0.f, FLOAT_BEFORE_ONE);
    Vec3f tangent, bitangent;
    BuildOrthonormalBasis(meanDirection, tangent, bitangent);
    return Normalize(tangent * (sinTheta * Cos(phi)) + bitangent * (sinTheta * Sin(phi)) +
                     meanDirection * cosTheta);
}

YBI_GPU Vec3f SampleVMM(const VMMDistribution &distribution, const Vec2f &u)
{
    YBI_ASSERT(distribution.numComponents > 0 &&
               distribution.numComponents <= VMM_MAX_NUM_COMPONENTS);
    const float random = Clamp(u.y, 0.f, FLOAT_BEFORE_ONE);
    float cdf = 0.f;
    float componentCdf = 0.f;
    uint32_t component = 0;
    for (uint32_t k = 0; k < distribution.numComponents; ++k)
    {
        const float weight = distribution.weights[k];
        if (weight > 0.f)
        {
            component = k;
            componentCdf = cdf;
            if (random < cdf + weight)
                break;
        }
        cdf += weight;
    }
    const float weight = distribution.weights[component];
    YBI_ASSERT(weight > 0.f);
    const float remapped = Clamp((random - componentCdf) / weight, 0.f, FLOAT_BEFORE_ONE);
    return SampleVMF(distribution.meanDirection[component],
                     distribution.kappa[component],
                     Vec2f(u.x, remapped));
}

// The PDF includes every lobe, including lobes that did not generate the sample.
YBI_GPU float SampleVMMPDF(const VMMDistribution &distribution, const Vec2f &u, Vec3f &direction)
{
    direction = SampleVMM(distribution, u);
    return EvalVMM(direction, distribution);
}

YBI_GPU bool SampleGuidingDirection(const PathGuidingKdTreeParams &params,
                                    const Vec3f &position,
                                    const Vec2f &u,
                                    Vec3f &direction,
                                    float &pdf)
{
    pdf = 0.f;
    VMMDistribution distribution;
    if (!PrepareGuidingDistribution(params, position, distribution))
        return false;
    pdf = SampleVMMPDF(distribution, u, direction);
    return pdf > 0.f && IsFinite(pdf);
}

YBI_GPU float
GuidingPDF(const PathGuidingKdTreeParams &params, const Vec3f &position, const Vec3f &direction)
{
    VMMDistribution distribution;
    return PrepareGuidingDistribution(params, position, distribution)
               ? EvalVMM(direction, distribution)
               : 0.f;
}

} // namespace ybi
