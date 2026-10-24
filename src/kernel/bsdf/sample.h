#pragma once

#include "core/vec3f.h"
#include "device/cuda_compat.h"

namespace ybi
{

enum BsdfSampleFlags : uint32_t
{
    BSDF_REFLECTION = 1u << 0,
    BSDF_TRANSMISSION = 1u << 1,
    BSDF_DELTA = 1u << 2,
    BSDF_DIFFUSE = 1u << 3,
    BSDF_GLOSSY = 1u << 4,
};

struct BsdfEvaluationProperties
{
    uint32_t flags = 0;
    float maxRoughness = -1.f;
};

struct BsdfSample
{
    Vec3f direction;
    Vec3f weight; // Cosine-weighted BSDF divided by the sampling probability.
    float pdf; // Solid-angle density, or the sampled delta branch's probability.
    float stochasticPdf;
    float eta; // Transmitted IOR / incident IOR; one for reflection.
    uint32_t flags;
    float sampledRoughness = -1.f; // Negative when roughness is unknown.
    float refractionDistanceScale = 1.f;
};

} // namespace ybi
