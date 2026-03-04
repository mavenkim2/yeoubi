#pragma once

#if defined(WITH_EMBREE)

#include "util/base.h"

YBI_NAMESPACE_BEGIN

struct CPUDevice;

void RegisterCPUDefaultKernels(CPUDevice *device);

YBI_NAMESPACE_END

#endif
