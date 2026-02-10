#include <algorithm>
#include <cstdint>
#include <cstring>

#if defined(_WIN32)
#include <windows.h>
#else
#include <sys/mman.h>
#include <unistd.h>
#endif

#include "util/memory_view.h"

#define ARENA_RESERVE_SIZE (64ll * 1024ll * 1024ll)
#define ARENA_COMMIT_SIZE (64ll * 1024ll)

YBI_NAMESPACE_BEGIN

class MemoryArena
{
private:
    uint8_t *base;
    size_t used;
    size_t committed;
    size_t reserveSize;
    size_t alignment;
    size_t commitSize;

public:
    explicit MemoryArena()
    {
        reserveSize = ARENA_RESERVE_SIZE;
        commitSize = ARENA_COMMIT_SIZE;

        assert(reserveSize % commitSize == 0 &&
               "Reserve size should be multiple of commit size.\n");
        alignment = 16;
#if defined(_WIN32)
        base = (uint8_t *)VirtualAlloc(nullptr, reserveSize, MEM_RESERVE, PAGE_NOACCESS);
#else
        m_base =
            (uint8_t *)mmap(nullptr, m_reserveSize, PROT_NONE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
#endif
        used = 0;
        committed = 0;
    }

    ~MemoryArena()
    {
        if (base)
        {
#if defined(_WIN32)
            VirtualFree(base, 0, MEM_RELEASE);
#else
            munmap(m_base, m_reserveSize);
#endif
        }
    }

    template <typename T>
    MemoryView<T> Push()
    {
        return PushArray(1);
    }

    template <typename T>
    MemoryView<T> PushArray(size_t count)
    {
        size_t bytesNeeded = count * sizeof(T);

        size_t alignedUsed = (used + alignment - 1) & ~(alignment - 1);
        size_t newUsed = alignedUsed + bytesNeeded;

        assert(newUsed <= reserveSize && "Arena out of address space!");

        if (newUsed > committed)
        {
            size_t commitTarget = (newUsed + commitSize - 1) & ~(commitSize - 1);
            size_t commitSize = commitTarget - committed;

#if defined(_WIN32)
            VirtualAlloc(base + committed, commitSize, MEM_COMMIT, PAGE_READWRITE);
#else
            mprotect(m_base + m_committed, commitSize, PROT_READ | PROT_WRITE);
#endif
            committed = commitTarget;
        }

        T *dest = (T *)(base + alignedUsed);

        used = newUsed;
        return {dest, count};
    }

    void Clear()
    {
        used = 0;
    }
};

YBI_NAMESPACE_END
