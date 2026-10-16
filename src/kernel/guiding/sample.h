#pragma once

#include "core/simd_float4.h"
#include "kernel/guiding/types.h"

namespace ybi
{

template <typename Float>
YBI_GPU Float LoadRadiance(const PathGuidingKdTreeParams &params,
                           uint32_t sampleIndex,
                           uint32_t numSamples);

// numValidSamples is a lane count, independent of the absolute sampleIndex.
template <typename Float>
YBI_GPU Float LoadRadianceContribution(const PathGuidingKdTreeParams &params,
                                       uint32_t sampleIndex,
                                       uint32_t numValidSamples);

template <typename Float>
YBI_GPU Vector3<Float>
LoadSampleDirection(const PathGuidingKdTreeParams &params, uint32_t sampleIndex, uint32_t numValidSamples);

template <>
YBI_GPU Vec3f LoadSampleDirection<float>(const PathGuidingKdTreeParams &params,
                                         uint32_t sampleIndex,
                                         uint32_t numValidSamples)
{
    YBI_ASSERT(numValidSamples == 1);
    return params.sampleDirection[sampleIndex];
}

template <>
YBI_GPU float LoadRadianceContribution<float>(const PathGuidingKdTreeParams &params,
                                              uint32_t sampleIndex,
                                              uint32_t numValidSamples)
{
    YBI_ASSERT(numValidSamples == 1);
    YBI_ASSERT(params.samplePdf[sampleIndex] > 0.f);
    return params.sampleIncidentRadianceEstimate[sampleIndex] / params.samplePdf[sampleIndex];
}

template <typename Sample>
YBI_GPU Sample LoadSampleData(const PathGuidingKdTreeParams &params,
                              uint32_t sampleIndex,
                              uint32_t numSamples);

template <>
YBI_GPU float LoadRadiance(const PathGuidingKdTreeParams &params,
                           uint32_t sampleIndex,
                           uint32_t numSamples)
{
    YBI_ASSERT(numSamples == 1);
    return params.sampleIncidentRadianceEstimate[sampleIndex];
}

template <>
YBI_GPU SampleData LoadSampleData<SampleData>(const PathGuidingKdTreeParams &params,
                                              uint32_t sampleIndex,
                                              uint32_t numSamples)
{
    YBI_ASSERT(numSamples == 1);
    SampleData sampleData{};
    sampleData.direction = params.sampleDirection[sampleIndex];
    sampleData.incidentRadianceEstimate = params.sampleIncidentRadianceEstimate[sampleIndex];
    sampleData.pdf = params.samplePdf[sampleIndex];
    return sampleData;
}

#if defined(SIMD_4)

template <>
inline Vector3<Float4> LoadSampleDirection<Float4>(const PathGuidingKdTreeParams &params,
                                                   uint32_t sampleIndex,
                                                   uint32_t numValidSamples)
{
    YBI_ASSERT(numValidSamples > 0 && numValidSamples <= Float4::Width);
    if (numValidSamples == Float4::Width)
        return LoadUnaligned(params.sampleDirection + sampleIndex);

    Vec3f directions[Float4::Width] = {};
    for (uint32_t lane = 0; lane < numValidSamples; ++lane)
        directions[lane] = params.sampleDirection[sampleIndex + lane];
    return LoadUnaligned(directions);
}

template <>
inline Float4 LoadRadianceContribution<Float4>(const PathGuidingKdTreeParams &params,
                                               uint32_t sampleIndex,
                                               uint32_t numValidSamples)
{
    YBI_ASSERT(numValidSamples > 0 && numValidSamples <= Float4::Width);
    Float4 radiance;
    Float4 pdf;
    if (numValidSamples == Float4::Width)
    {
        radiance = LoadUnaligned(params.sampleIncidentRadianceEstimate + sampleIndex);
        pdf = LoadUnaligned(params.samplePdf + sampleIndex);
    }
    else
    {
        // A short range need not have padding; never read beyond its end.
        float radianceLanes[Float4::Width] = {};
        float pdfLanes[Float4::Width] = {1.f, 1.f, 1.f, 1.f};
        for (uint32_t lane = 0; lane < numValidSamples; ++lane)
        {
            radianceLanes[lane] = params.sampleIncidentRadianceEstimate[sampleIndex + lane];
            pdfLanes[lane] = params.samplePdf[sampleIndex + lane];
        }
        radiance = LoadUnaligned(radianceLanes);
        pdf = LoadUnaligned(pdfLanes);
    }
    YBI_ASSERT(All(pdf > Float4(0.f)));
    return radiance / pdf;
}

template <>
YBI_GPU Float4 LoadRadiance(const PathGuidingKdTreeParams &params,
                            uint32_t sampleIndex,
                            uint32_t numSamples)
{
    YBI_ASSERT(numSamples > 0 && numSamples <= Float4::Width);
    if (numSamples == Float4::Width)
        return LoadUnaligned(params.sampleIncidentRadianceEstimate + sampleIndex);

    float radiance[Float4::Width] = {};
    for (uint32_t lane = 0; lane < numSamples; ++lane)
        radiance[lane] = params.sampleIncidentRadianceEstimate[sampleIndex + lane];
    return LoadUnaligned(radiance);
}

template <>
inline VMMSampleData<Float4> LoadSampleData<VMMSampleData<Float4>>(const PathGuidingKdTreeParams &params,
                                                                   uint32_t sampleIndex,
                                                                   uint32_t numSamples)
{
    YBI_ASSERT(numSamples > 0 && numSamples <= Float4::Width);
    VMMSampleData<Float4> sample{};
    sample.direction = LoadSampleDirection<Float4>(params, sampleIndex, numSamples);
    sample.incidentRadianceEstimate = LoadRadiance<Float4>(params, sampleIndex, numSamples);
    if (numSamples == Float4::Width)
    {
        sample.pdf = LoadUnaligned(params.samplePdf + sampleIndex);
    }
    else
    {
        float pdf[Float4::Width] = {1.f, 1.f, 1.f, 1.f};
        for (uint32_t lane = 0; lane < numSamples; ++lane)
            pdf[lane] = params.samplePdf[sampleIndex + lane];
        sample.pdf = LoadUnaligned(pdf);
    }
    return sample;
}
#endif

} // namespace ybi
