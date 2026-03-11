#pragma once

#if defined(_WIN32)
#include <windows.h>
#else
#include <sys/mman.h>
#include <unistd.h>
#endif
#include "util/memory_arena.h"
#include "util/memory_view.h"

namespace ybi
{

struct CPUMemoryProvider
{
    template <typename T>
    using View = MemoryView<T>;

    static uint8_t *Reserve(size_t size)
    {
#if defined(_WIN32)
        return (uint8_t *)VirtualAlloc(nullptr, size, MEM_RESERVE, PAGE_NOACCESS);
#else
        return (uint8_t *)mmap(nullptr, size, PROT_NONE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
#endif
    }
    static void *Commit(uint8_t *ptr, size_t size)
    {
#if defined(_WIN32)
        VirtualAlloc(ptr, size, MEM_COMMIT, PAGE_READWRITE);
#else
        mprotect(ptr, size, PROT_READ | PROT_WRITE);
#endif
        return nullptr;
    }
    static void
    Free(uint8_t *base, size_t reserveSize, const ArenaCommit *commits, size_t numCommits)
    {
        (void)commits;
        (void)numCommits;
#if defined(_WIN32)
        VirtualFree(base, 0, MEM_RELEASE);
#else
        munmap(base, reserveSize);
#endif
    }
    static size_t GetPageSize()
    {
#if defined(_WIN32)
        SYSTEM_INFO si;
        GetSystemInfo(&si);
        return (size_t)si.dwPageSize;
#else
        return (size_t)sysconf(_SC_PAGESIZE);
#endif
    }
};

using HostMemoryArena = MemoryArena<CPUMemoryProvider>;

} // namespace ybi
