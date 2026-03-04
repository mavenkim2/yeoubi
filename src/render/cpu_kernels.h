#pragma once

#if defined(WITH_EMBREE)

#include "render/dispatch_types.h"
#include "util/base.h"

YBI_NAMESPACE_BEGIN

bool CPUDispatchKernel(const DispatchParams &dispatchParams, RenderKernelId kernelId);

YBI_NAMESPACE_END

#endif
