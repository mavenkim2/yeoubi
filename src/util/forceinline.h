#pragma once

#if !defined(_MSC_VER) && !defined(__forceinline)
#if defined(__GNUC__) || defined(__clang__)
#define __forceinline inline __attribute__((always_inline))
#else
#define __forceinline inline
#endif
#endif
