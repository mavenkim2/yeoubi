#pragma once

#include "core/math.h"
#include "kernel/bsdf/sample.h"
#include "kernel/guiding/sampling.h"

namespace ybi
{

template <typename BsdfSampleFunc, typename BsdfEvalFunc>
YBI_GPU bool SampleGuiding(BsdfSample &sample,
                           const PathGuidingKdTreeParams *params,
                           const Vec3f &position,
                           const BsdfSampleFunc &sampleBsdf,
                           const BsdfEvalFunc &evalBsdf,
                           const Vec3f &u)
{
    VMMDistribution distribution;
    if (!params || !PrepareGuidingDistribution(*params, position, distribution))
    {
        const bool success = sampleBsdf(sample);
        if (success)
            sample.stochasticPdf = sample.pdf;
        return success;
    }

    constexpr float alpha = 0.5f;
    constexpr float defensiveConstant = 1.f / (4.f * YBI_PI);

    struct RISSample
    {
        float candidatePdf = 0.f;
        float targetFunction = 0.f;
        float candidateResamplingWeight = 0.f;

        YBI_GPU RISSample(const BsdfSample &bsdf, bool success, float guidingPdf)
        {
            if (success && bsdf.weight.MaxComponent() != 0.f && bsdf.pdf != 0.f)
            {
                const float bsdfCos = (bsdf.weight * bsdf.pdf).MaxComponent();
                candidatePdf = alpha * bsdf.pdf + (1.f - alpha) * guidingPdf;
                targetFunction =
                    bsdfCos * (alpha * defensiveConstant + (1.f - alpha) * guidingPdf);
                candidateResamplingWeight = targetFunction / candidatePdf;
            }
        }
    };

    BsdfSample bsdfSample{};
    const bool success = sampleBsdf(bsdfSample);

    Vec3f sampledGuidingDirection;
    const float sampledGuidingPdf = SampleVMMPDF(distribution, u.Xy(), sampledGuidingDirection);
    BsdfSample guidingSample = evalBsdf(sampledGuidingDirection);

    const float bsdfGuidingPdf = success && !(bsdfSample.flags & BSDF_DELTA)
                                     ? EvalVMM(bsdfSample.direction, distribution)
                                     : 0.f;
    RISSample bsdfRisSample(bsdfSample, success, bsdfGuidingPdf);
    RISSample guidingRisSample(guidingSample, true, sampledGuidingPdf);

    const float bsdfWeight = bsdfRisSample.candidateResamplingWeight;
    const float guidingWeight = guidingRisSample.candidateResamplingWeight;
    const float totalWeight = bsdfWeight + guidingWeight;

    if (totalWeight == 0.f)
    {
        return false;
    }

    const float rand = u.z;
    const float bsdfResampleProb = bsdfWeight / totalWeight;
    const float guidingResampleProb = guidingWeight / totalWeight;

    const bool selectedBsdf =
        guidingWeight == 0.f || (bsdfWeight > 0.f && rand * totalWeight < bsdfWeight);
    if (selectedBsdf)
    {
        const Vec3f bsdfCos = bsdfSample.weight * bsdfSample.pdf;
        sample = bsdfSample;
        sample.direction = bsdfSample.direction;
        sample.stochasticPdf = 2.f * bsdfRisSample.candidatePdf * bsdfResampleProb;
        sample.weight = bsdfCos / sample.stochasticPdf;
        sample.pdf = bsdfRisSample.candidatePdf;
    }
    else
    {
        const Vec3f guidingBsdfCos = guidingSample.weight * guidingSample.pdf;
        sample = guidingSample;
        sample.direction = sampledGuidingDirection;
        sample.stochasticPdf = 2.f * guidingRisSample.candidatePdf * guidingResampleProb;
        sample.weight = guidingBsdfCos / sample.stochasticPdf;
        sample.pdf = guidingRisSample.candidatePdf;
    }

    return true;
}

// The same 50/50 proposal PDF is used by light MIS.
YBI_GPU float GuidingProposalPDF(const PathGuidingKdTreeParams *params,
                                 const Vec3f &position,
                                 const Vec3f &direction,
                                 float bsdfPdf)
{
    VMMDistribution distribution;
    if (!params || !PrepareGuidingDistribution(*params, position, distribution))
        return bsdfPdf;
    return 0.5f * bsdfPdf + 0.5f * EvalVMM(direction, distribution);
}

} // namespace ybi
