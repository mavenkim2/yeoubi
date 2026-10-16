#pragma once

#include "device/cuda_compat.h"

namespace ybi
{

template <typename T, uint32_t BlockSize = 256>
T BlockReduceSum(T *values)
{
    static_assert(BlockSize >= 32 && BlockSize <= 1024 && BlockSize % 32 == 0,
                  "Block reduction requires complete warps");
    T warpValues[32] = {};
    for (uint32_t warp = 0; warp < BlockSize / 32; ++warp)
    {
        T *lanes = values + warp * 32;
        for (uint32_t offset = 16; offset != 0; offset >>= 1)
            for (uint32_t lane = 0; lane < offset; ++lane)
                lanes[lane] = lanes[lane] + lanes[lane + offset];
        warpValues[warp] = lanes[0];
    }
    for (uint32_t offset = 16; offset != 0; offset >>= 1)
        for (uint32_t lane = 0; lane < offset; ++lane)
            warpValues[lane] = warpValues[lane] + warpValues[lane + offset];
    return warpValues[0];
}

} // namespace ybi
