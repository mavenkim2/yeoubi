#pragma once

#include "kernel/guiding/guiding.h"

namespace ybi
{

struct VMMUpdateParams
{
    PathGuidingKdTreeParams samples;
    VMMDistribution *distributions;
    // Initialized, normalized previous statistics; updated once per batch.
    SufficientStatistics *statistics;
    const VMMUpdateRange *ranges;
    uint32_t numDistributions;
    VMMFitConfiguration configuration;
};

} // namespace ybi
