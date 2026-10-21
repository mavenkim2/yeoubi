#pragma once

#include "device/cuda_compat.h"

#define YBI_KERNEL extern "C" __global__ void
#define YBI_ITERATE_WORK_ITEMS(index, count) \
    for (uint32_t index = blockIdx.x; index < (count); index += gridDim.x)
#define YBI_ITERATE_ELEMENTS(index, count) \
    for (uint32_t index = blockIdx.x * blockDim.x + threadIdx.x; \
         index < (count); index += blockDim.x * gridDim.x)
#define YBI_END_ITERATE_WORK_ITEMS
#define YBI_FIRST_THREAD (blockIdx.x == 0 && threadIdx.x == 0)
