#pragma once

#include "util/assert.h"

YBI_NAMESPACE_BEGIN

#ifdef WITH_CUDA
#ifdef YBI_DEBUG
#define CUDA_ASSERT(statement) \
    { \
        CUresult result = statement; \
        if (result != CUDA_SUCCESS) \
        { \
            const char *str, *name; \
            cuGetErrorString(result, &str); \
            cuGetErrorName(result, &name); \
            fprintf(stderr, \
                    "CUDA Error (%s): %s in %s (%s:%d)\n", \
                    name, \
                    str, \
                    #statement, \
                    __FILE__, \
                    __LINE__); \
            BREAK_IN_DEBUGGER(); \
        } \
    }

#else
#define CUDA_ASSERT(statement) statement
#endif

#ifdef WITH_OPTIX
#ifdef YBI_DEBUG

#define OPTIX_ASSERT(statement) \
    { \
        OptixResult result = statement; \
        if (result != OPTIX_SUCCESS) \
        { \
            const char *str = optixGetErrorString(result); \
            const char *name = optixGetErrorName(result); \
            fprintf(stderr, \
                    "Optix Error (%s): %s in %s (%s:%d)\n", \
                    name, \
                    str, \
                    #statement, \
                    __FILE__, \
                    __LINE__); \
            BREAK_IN_DEBUGGER(); \
        } \
    }
#else
#define OPTIX_ASSERT(statement) statement
#endif

#endif
#endif

YBI_NAMESPACE_END
