#pragma once

#include "core/octahedral.h"
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

template <typename Float>
YBI_GPU Float LoadSampleDistance(const PathGuidingKdTreeParams &params,
                                 uint32_t sampleIndex,
                                 uint32_t numValidSamples);

template <>
YBI_GPU float LoadSampleDistance<float>(const PathGuidingKdTreeParams &params,
                                        uint32_t sampleIndex,
                                        uint32_t numValidSamples)
{
    YBI_ASSERT(numValidSamples == 1);
    return params.sampleDistance ? params.sampleDistance[sampleIndex] : FLOAT_INFINITY;
}

template <>
YBI_GPU Vec3f LoadSampleDirection<float>(const PathGuidingKdTreeParams &params,
                                         uint32_t sampleIndex,
                                         uint32_t numValidSamples)
{
    YBI_ASSERT(numValidSamples == 1);
    return DecodeOctahedral(params.sampleDirection[sampleIndex]);
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
    sampleData.direction = LoadSampleDirection<float>(params, sampleIndex, numSamples);
    sampleData.incidentRadianceEstimate = params.sampleIncidentRadianceEstimate[sampleIndex];
    sampleData.pdf = params.samplePdf[sampleIndex];
    sampleData.distance = LoadSampleDistance<float>(params, sampleIndex, numSamples);
    return sampleData;
}

#if defined(SIMD_4)

template <>
inline Float4 LoadSampleDistance<Float4>(const PathGuidingKdTreeParams &params,
                                         uint32_t sampleIndex,
                                         uint32_t numValidSamples)
{
    YBI_ASSERT(numValidSamples > 0 && numValidSamples <= Float4::Width);
    if (numValidSamples == Float4::Width)
        return params.sampleDistance ? LoadUnaligned(params.sampleDistance + sampleIndex)
                                     : Float4(FLOAT_INFINITY);

    float distances[Float4::Width] = {1.f, 1.f, 1.f, 1.f};
    for (uint32_t lane = 0; lane < numValidSamples; ++lane)
        distances[lane] = params.sampleDistance ? params.sampleDistance[sampleIndex + lane]
                                               : FLOAT_INFINITY;
    return LoadUnaligned(distances);
}

template <>
inline Vector3<Float4> LoadSampleDirection<Float4>(const PathGuidingKdTreeParams &params,
                                                   uint32_t sampleIndex,
                                                   uint32_t numValidSamples)
{
    YBI_ASSERT(numValidSamples > 0 && numValidSamples <= Float4::Width);
    if (numValidSamples == Float4::Width)
        return DecodeOctahedral4(params.sampleDirection + sampleIndex);

    // A short range need not have padding; decode only after copying its valid entries.
    uint32_t directions[Float4::Width] = {};
    for (uint32_t lane = 0; lane < numValidSamples; ++lane)
        directions[lane] = params.sampleDirection[sampleIndex + lane];
    const auto decoded = DecodeOctahedral4(directions);
    const Mask4 valid = Mask4::Create((1u << numValidSamples) - 1u);
    return {Select(valid, decoded.x, Float4(0.f)),
            Select(valid, decoded.y, Float4(0.f)),
            Select(valid, decoded.z, Float4(0.f))};
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
    sample.distance = LoadSampleDistance<Float4>(params, sampleIndex, numSamples);
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
