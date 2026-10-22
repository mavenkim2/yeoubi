#include "kernel/backend/cuda/kernel.h"
#include "kernel/backend/cuda/guiding/vmm_backend.h"

#define YBI_GUIDING_BACKEND CudaVMMBackend
#include "kernel/guiding/dispatch.h"
#undef YBI_GUIDING_BACKEND
#undef YBI_KERNEL
#undef YBI_SHARED
#undef YBI_ITERATE_WORK_ITEMS
#undef YBI_ITERATE_ELEMENTS
#undef YBI_END_ITERATE_WORK_ITEMS
#undef YBI_FIRST_THREAD
