#pragma once

#include <cstdint>
#include <vector>

namespace ybi
{

float CreatePiecewise1DCDF(const float *values, uint32_t numValues, float *cdf);
float CreatePiecewise2DCDF(const std::vector<float> &values,
                          uint32_t width,
                          uint32_t height,
                          std::vector<float> &conditionalCDF,
                          std::vector<float> &marginalCDF);

} // namespace ybi
