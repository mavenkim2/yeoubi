#pragma once

#include <algorithm>
#include <cstdint>
#include <cstring>

#if defined(_WIN32)
#include <windows.h>
#else
#include <sys/mman.h>
#include <unistd.h>
#endif

#include "util/aligned_malloc.h"
#include "util/base.h"
#include "util/memory_view.h"

#define ARENA_RESERVE_SIZE (64ll * 1024ll * 1024ll)
#define ARENA_COMMIT_SIZE (64ll * 1024ll)

YBI_NAMESPACE_BEGIN

inline size_t GetPageSize()
{
#if defined(_WIN32)
    SYSTEM_INFO si;
    GetSystemInfo(&si);
    return (size_t)si.dwPageSize;
#else
    return (size_t)sysconf(_SC_PAGESIZE);
#endif
}

class MemoryArena
{
private:
    MemoryArena *current;
    MemoryArena *prev;
    uint8_t *base;
    size_t used;
    size_t committed;
    size_t reserveSize;
    size_t alignment;
    size_t commitSize;

public:
    explicit MemoryArena(size_t requestedReserveSize = ARENA_RESERVE_SIZE)
    {
        size_t pageSize = GetPageSize();
        reserveSize = util::AlignUp(requestedReserveSize, pageSize);
        commitSize = util::AlignUp(ARENA_COMMIT_SIZE, pageSize);

        assert(reserveSize % commitSize == 0 &&
               "Reserve size should be multiple of commit size.\n");
        alignment = 16;
#if defined(_WIN32)
        base = (uint8_t *)VirtualAlloc(nullptr, reserveSize, MEM_RESERVE, PAGE_NOACCESS);
#else
        base =
            (uint8_t *)mmap(nullptr, reserveSize, PROT_NONE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
#endif
        used = 0;
        committed = 0;
        prev = nullptr;
        current = this;

        assert(base != nullptr && "Virtual Memory Reservation Failed");
    }

    ~MemoryArena()
    {
        if (current && current != this)
        {
            MemoryArena *curr = current;
            while (curr != this)
            {
                MemoryArena *toDelete = curr;
                curr = curr->prev;
                toDelete->prev = nullptr;
                toDelete->~MemoryArena();
                util::AlignedFree(toDelete);
            }
        }
        if (base)
        {
#if defined(_WIN32)
            VirtualFree(base, 0, MEM_RELEASE);
#else
            munmap(base, reserveSize);
#endif
        }
    }

    MemoryArena(const MemoryArena &) = delete;
    MemoryArena &operator=(const MemoryArena &) = delete;
    MemoryArena(MemoryArena &&) = delete;
    MemoryArena &operator=(MemoryArena &&) = delete;

    template <typename T>
    MemoryView<T> Push()
    {
        return PushArray(1);
    }

    template <typename T>
    MemoryView<T> PushArray(size_t count)
    {
        size_t bytesNeeded = count * sizeof(T);

        size_t alignedUsed = util::AlignUp(current->used, current->alignment);
        size_t newUsed = alignedUsed + bytesNeeded;

        if (newUsed > current->reserveSize)
        {
            size_t newBlockSize = current->reserveSize;
            if (bytesNeeded > current->reserveSize)
            {
                newBlockSize = bytesNeeded;
            }

            void *newPtr = util::AlignedAlloc(sizeof(MemoryArena), 8);
            MemoryArena *newArena = new (newPtr) MemoryArena(newBlockSize);
            newArena->prev = current;
            current = newArena;

            alignedUsed = 0;
            newUsed = bytesNeeded;
        }

        if (newUsed > current->committed)
        {
            size_t commitTarget = util::AlignUp(newUsed, current->commitSize);
            size_t commitRequestSize = commitTarget - current->committed;

#if defined(_WIN32)
            VirtualAlloc(
                current->base + current->committed, commitRequestSize, MEM_COMMIT, PAGE_READWRITE);
#else
            mprotect(
                current->base + current->committed, commitRequestSize, PROT_READ | PROT_WRITE);
#endif
            current->committed = commitTarget;
        }

        T *dest = (T *)(current->base + alignedUsed);

        current->used = newUsed;
        return {dest, count};
    }
    void Clear()
    {
        MemoryArena *curr = current;
        while (curr != this)
        {
            MemoryArena *toDelete = curr;
            curr = curr->prev;
            toDelete->prev = nullptr;
            toDelete->~MemoryArena();
            util::AlignedFree(toDelete);
            delete toDelete;
        }
        current = this;
        used = 0;
    }
};

YBI_NAMESPACE_END
