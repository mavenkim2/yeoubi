#pragma once

#include "core/octahedral.h"
#include "kernel/guiding/types.h"

namespace ybi
{

YBI_GPU const PathGuidingTrainingSample &LoadTrainingSample(
    const PathGuidingKdTreeParams &params, uint32_t sampleIndex)
{
    YBI_ASSERT(sampleIndex < params.numSampleIndices);
    return params.sampleData[params.sampleIndices[sampleIndex]];
}

template <typename Float>
YBI_GPU Float LoadRadiance(const PathGuidingKdTreeParams &params,
                           uint32_t sampleIndex,
                           uint32_t numSamples);

template <typename Float>
YBI_GPU Float LoadRadianceContribution(const PathGuidingKdTreeParams &params,
                                       uint32_t sampleIndex,
                                       uint32_t numValidSamples);

template <typename Float>
YBI_GPU Vector3<Float> LoadSampleDirection(const PathGuidingKdTreeParams &params,
                                           uint32_t sampleIndex,
                                           uint32_t numValidSamples);

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
    return LoadTrainingSample(params, sampleIndex).distance;
}

template <>
YBI_GPU Vec3f LoadSampleDirection<float>(const PathGuidingKdTreeParams &params,
                                         uint32_t sampleIndex,
                                         uint32_t numValidSamples)
{
    YBI_ASSERT(numValidSamples == 1);
    return DecodeOctahedral(LoadTrainingSample(params, sampleIndex).octahedralDirection);
}

template <>
YBI_GPU float LoadRadianceContribution<float>(const PathGuidingKdTreeParams &params,
                                              uint32_t sampleIndex,
                                              uint32_t numValidSamples)
{
    YBI_ASSERT(numValidSamples == 1);
    const PathGuidingTrainingSample &sample = LoadTrainingSample(params, sampleIndex);
    YBI_ASSERT(sample.pdf > 0.f);
    return sample.GetIncidentRadiance() / sample.pdf;
}

template <>
YBI_GPU float LoadRadiance<float>(const PathGuidingKdTreeParams &params,
                                  uint32_t sampleIndex,
                                  uint32_t numSamples)
{
    YBI_ASSERT(numSamples == 1);
    return LoadTrainingSample(params, sampleIndex).GetIncidentRadiance();
}

template <typename Sample>
YBI_GPU Sample LoadSampleData(const PathGuidingKdTreeParams &params,
                              uint32_t sampleIndex,
                              uint32_t numSamples);

template <>
YBI_GPU SampleData LoadSampleData<SampleData>(const PathGuidingKdTreeParams &params,
                                              uint32_t sampleIndex,
                                              uint32_t numSamples)
{
    YBI_ASSERT(numSamples == 1);
    const PathGuidingTrainingSample &sample = LoadTrainingSample(params, sampleIndex);
    SampleData result{};
    result.direction = DecodeOctahedral(sample.octahedralDirection);
    result.position = sample.position;
    result.incidentRadianceEstimate = sample.GetIncidentRadiance();
    result.pdf = sample.pdf;
    result.distance = sample.distance;
    return result;
}

} // namespace ybi
