#pragma once

#include "kernel/sampling/types.h"
#include <vector>

namespace ybi
{

float CreateAliasTable(const float *weights, uint32_t count, AliasEntry *table);
float CreateAliasTable2D(const std::vector<float> &weights,
                         uint32_t width,
                         uint32_t height,
                         std::vector<AliasEntry> &conditionalAlias,
                         std::vector<AliasEntry> &marginalAlias);

} // namespace ybi
