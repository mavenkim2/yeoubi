// Implementation adapted from https://github.com/OpenPathGuidingLibrary/openpgl.
// Copyright 2020 Intel Corporation
// SPDX-License-Identifier: Apache-2.0
// Full license: licenses/OpenPGL.txt (relative to the repository root).
#pragma once

#include "core/lane.h"
#include "core/math.h"
#include "core/vec3f.h"
#include "device/cuda_compat.h"
#include "kernel/guiding/types.h"
#include <cmath>

namespace ybi
{

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

YBI_CPU_GPU float VMFNormalization(float kappa)
{
    return kappa == 0.f ? 1.f / (4.f * YBI_PI) : kappa / (2.f * YBI_PI * -expm1f(-2.f * kappa));
}

// TODO: maybe validation function that makes sure weights sum to 1

template <typename Float>
YBI_GPU Float EvalVMF(const Vector3<Float> &w, const Vec3f &meanDirection, float kappa)
{
    if (kappa == 0.f)
    {
        return Float(1.f / (4.f * YBI_PI));
    }

    const Float cosine =
        w.x * Float(meanDirection.x) + w.y * Float(meanDirection.y) + w.z * Float(meanDirection.z);
    const Float exponent = Float(kappa) * (cosine - Float(1.f));
    return Float(VMFNormalization(kappa)) * Exp(exponent);
}

template <typename Float>
YBI_GPU bool SoftAssignment(const VMMDistribution &distribution,
                            const Vector3<Float> &sampleDirection,
                            Float *softAssignment,
                            Float &pdf)
{
    Float V(0.f);

    for (uint32_t componentIndex = 0; componentIndex < distribution.numComponents;
         componentIndex++)
    {
        const Float v = Float(distribution.weights[componentIndex]) *
                        EvalVMF<Float>(sampleDirection,
                                       distribution.meanDirection[componentIndex],
                                       distribution.kappa[componentIndex]);
        V += v;

        softAssignment[componentIndex] = v;
    }

    pdf = V;
    if (pdf <= Float(1e-16f))
        return false;
    for (uint32_t k = 0; k < distribution.numComponents; ++k)
    {
        softAssignment[k] /= pdf;
    }
    return true;
}

YBI_GPU float EvalVMF(const Vec3f &w, const Vec3f &meanDirection, float kappa)
{
    return EvalVMF<float>(w, meanDirection, kappa);
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
