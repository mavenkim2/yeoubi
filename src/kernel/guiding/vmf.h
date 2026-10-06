#pragma once

#include "core/math.h"
#include "core/vec3f.h"
#include "device/cuda_compat.h"
#include <cmath>

namespace ybi
{

#define VMM_MAX_NUM_COMPONENTS 16

struct VMMDistribution
{
    Vec3f meanDirection[VMM_MAX_NUM_COMPONENTS];
    // NOTE: must sum to 1
    float weights[VMM_MAX_NUM_COMPONENTS];
    float kappa[VMM_MAX_NUM_COMPONENTS];

    uint32_t numComponents;
};

YBI_CPU_GPU float KappaToMeanCosine(float kappa)
{
    if (kappa <= 0.f)
    {
        return 0.f;
    }
    return 1.f / tanhf(kappa) - 1.f / kappa;
}

YBI_CPU_GPU float MeanCosineToKappa(const float &meanCosine)
{
    const float ones = 1.f;
    const float dim = 3.f;
    const float meanCosine2 = meanCosine * meanCosine;
    return (meanCosine * dim - meanCosine * meanCosine2) / (ones - meanCosine2);
}

// TODO: maybe validation function that makes sure weights sum to 1

YBI_GPU float EvalVMF(const Vec3f &w, const Vec3f &meanDirection, float kappa)
{
    if (kappa == 0.f)
    {
        return 1.f / (4.f * YBI_PI);
    }

    float exponent = kappa * (Dot(w, meanDirection) - 1.f);
    float divisor = 2.f * YBI_PI * (1.f - Exp(-2.f * kappa));
    float result = kappa / divisor * Exp(exponent);
    return result;
}

YBI_GPU float EvalVMM(const Vec3f &w, const VMMDistribution &distribution)
{
    float result = 0.f;
    for (int componentIndex = 0; componentIndex < distribution.numComponents; componentIndex++)
    {
        result += distribution.weights[componentIndex] *
                  EvalVMF(w,
                          distribution.meanDirection[componentIndex],
                          distribution.kappa[componentIndex]);
    }
    return result;
}

} // namespace ybi
