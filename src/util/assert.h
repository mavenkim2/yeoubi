#pragma once

#include <cstdio>

#ifdef YBI_DEBUG
#if defined(_MSC_VER)
#include <intrin.h>
#define BREAK_IN_DEBUGGER() __debugbreak()
#elif defined(__GNUC__) || defined(__clang__)
#define BREAK_IN_DEBUGGER() __builtin_trap()
#else
#include <cstdlib>
#define BREAK_IN_DEBUGGER() std::abort()
#endif

#define LOG_ERROR_INTERNAL(type, expr, message) \
    fprintf(stderr, "[%s] %s:%d\nEXPR: %s\n%s", type, __FILE__, __LINE__, #expr, message)

#define ASSERT(expr) \
    do \
    { \
        if (!(expr)) \
        { \
            LOG_ERROR_INTERNAL("ASSERTION FAILED", expr, ""); \
            BREAK_IN_DEBUGGER(); \
        } \
    } while (0)

#define CHECK(expr, msg) \
    do \
    { \
        if (!(expr)) \
        { \
            LOG_ERROR_INTERNAL("CHECK FAILED", expr, msg); \
            BREAK_IN_DEBUGGER(); \
        } \
    } while (0)

#define ERROR(expr, msg) \
    do \
    { \
        if (!(expr)) \
        { \
            LOG_ERROR_INTERNAL("FATAL ERROR", expr, msg); \
            BREAK_IN_DEBUGGER(); \
        } \
    } while (0)

#define LOGFATAL(msg) \
    do \
    { \
        fprintf(stderr, "[%s] %s:%d\n%s", "FATAL ERROR", __FILE__, __LINE__, #msg); \
        BREAK_IN_DEBUGGER(); \
    } while (0)

#else
#define ASSERT(expr) ((void)(0))
#define CHECK(expr, msg) ((void)expr)
#define ERROR(expr, msg) ((void)(0))
#endif
