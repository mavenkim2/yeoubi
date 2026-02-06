#pragma once

YBI_NAMESPACE_BEGIN

namespace util
{

inline void *AlignedAlloc(size_t size, size_t alignment)
{
#if defined(_WIN32)
    return _aligned_malloc(size, alignment);
#else
    void *ptr = nullptr;
    if (posix_memalign(&ptr, alignment, size) != 0) return nullptr;
    return ptr;
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

} // namespace util

YBI_NAMESPACE_END
