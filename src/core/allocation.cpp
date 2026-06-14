#include "core/allocation.h"
#include "core/assert.h"
#include "core/bit.h"
#include <cstdlib>

#if defined(_WIN32)
#include <malloc.h>
#endif

namespace ybi
{

void *AlignedAlloc(size_t size, size_t alignment)
{
    YBI_ASSERT(size != 0 && IsPowerOfTwo(alignment));

#if defined(_WIN32)
    return _aligned_malloc(size, alignment);
#else
    void *ptr = nullptr;
    return posix_memalign(&ptr, alignment, size) == 0 ? ptr : nullptr;
#endif
}

void AlignedFree(void *ptr)
{
#if defined(_WIN32)
    _aligned_free(ptr);
#else
    free(ptr);
#endif
}

} // namespace ybi
