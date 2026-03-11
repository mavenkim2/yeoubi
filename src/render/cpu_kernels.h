#pragma once

#if defined(WITH_EMBREE)

#include "render/dispatch_types.h"
#include "util/base.h"

namespace ybi
{

bool CPUDispatchKernel(const DispatchParams &dispatchParams, RenderKernelId kernelId);

} // namespace ybi

#endif
