#pragma once

#include "core/assert.h"
#include "device/cuda_compat.h"
#include <tbb/parallel_for.h>

#define YBI_KERNEL void
#define YBI_SHARED
#define YBI_ITERATE_WORK_ITEMS(index, count) \
    tbb::parallel_for(0u, (count), [&](uint32_t index)
#define YBI_ITERATE_ELEMENTS(index, count) YBI_ITERATE_WORK_ITEMS(index, count)
#define YBI_END_ITERATE_WORK_ITEMS )
#define YBI_FIRST_THREAD true

namespace ybi
{

template <typename ParamsType, void (*Kernel)(ParamsType *)>
void DispatchKernelCPU(void **arguments)
{
    YBI_ASSERT(arguments && arguments[0]);
    ParamsType *params = *static_cast<ParamsType **>(arguments[0]);
    YBI_ASSERT(params);
    Kernel(params);
}

} // namespace ybi
