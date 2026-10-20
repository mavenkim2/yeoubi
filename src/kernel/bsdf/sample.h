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
};

struct BsdfSample
{
    Vec3f direction;
    Vec3f weight; // Cosine-weighted BSDF divided by the sampling probability.
    float pdf; // Solid-angle density, or the sampled delta branch's probability.
    float eta; // Transmitted IOR / incident IOR; one for reflection.
    uint32_t flags;
};

} // namespace ybi
