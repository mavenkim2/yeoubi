#pragma once

#include "device/cuda_compat.h"

namespace ybi
{

struct AliasEntry
{
    float threshold;
    uint32_t alias;
    float probability;
};

} // namespace ybi
