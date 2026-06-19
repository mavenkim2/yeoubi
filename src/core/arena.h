#pragma once

#include "core/allocation.h"
#include "core/array.h"
#include "core/assert.h"
#include "core/bit.h"
#include <cstdint>

#if defined(_WIN32)
#include <windows.h>
#else
#include <sys/mman.h>
#include <unistd.h>
#endif

#define ARENA_RESERVE_SIZE (64ll * 1024ll * 1024ll)
#define ARENA_COMMIT_SIZE (64ll * 1024ll)

namespace ybi
{

template <typename T>
struct View
{
    using value_type = T;

    T *ptr = nullptr;
    size_t count = 0;

    T &operator[](size_t index)
    {
        YBI_ASSERT(index < count);
        return ptr[index];
    }
    T &operator[](size_t index) const
    {
        YBI_ASSERT(index < count);
        return ptr[index];
    }

    T *begin() const
    {
        return ptr;
    }
    T *end() const
    {
        return ptr + count;
    }
    T *data() const
    {
        return ptr;
    }
    T *data()
    {
        return ptr;
    }
    size_t size() const
    {
        return count;
    }
    size_t numBytes() const
    {
        return count * sizeof(T);
    }

    View<T> &operator+=(size_t offset)
    {
        YBI_ASSERT(offset <= count);
        ptr += offset;
        count -= offset;
        return *this;
    }

    View<uint8_t> CastToBytes() const
    {
        View<uint8_t> view;
        view.ptr = (uint8_t *)ptr;
        view.count = numBytes();
        return view;
    }
};

class Arena
{
    struct ArenaCommit
    {
        uint8_t *ptr;
        size_t size;
        void *cookie;

        ArenaCommit(uint8_t *ptr, size_t size, void *cookie) : ptr(ptr), size(size), cookie(cookie)
        {
        }
    };

public:
    struct Mark
    {
        Arena *arena;
        Arena *current;
        size_t used;
    };

    struct Scope
    {
        Mark mark;

        explicit Scope(Arena *arena)
        {
            YBI_ASSERT(arena);
            mark.arena = arena;
            mark.current = arena->current;
            mark.used = arena->current->used;
        }
        ~Scope()
        {
            if (mark.arena)
            {
                Arena *curr = mark.arena->current;

                while (curr != mark.current)
                {
                    YBI_ASSERT(curr);
                    Arena *toDelete = curr;
                    curr = curr->prev;
                    toDelete->prev = nullptr;
                    toDelete->~Arena();
                    AlignedFree(toDelete);
                }
                mark.arena->current = mark.current;
                mark.current->used = mark.used;
            }
        }
        Scope(const Scope &) = delete;
        Scope &operator=(const Scope &) = delete;
    };

private:
    Arena *current;
    Arena *prev;
    uint8_t *base;
    size_t used;
    size_t committed;
    size_t reserveSize;
    size_t alignment;
    size_t commitSize;
    Array<ArenaCommit> commits;

public:
    explicit Arena(size_t requestedReserveSize = ARENA_RESERVE_SIZE)
    {
        size_t pageSize = GetPageSize();
        reserveSize = AlignUp(requestedReserveSize, pageSize);
        commitSize = AlignUp(ARENA_COMMIT_SIZE, pageSize);

        alignment = 16;

        base = Reserve(reserveSize);
        used = 0;
        committed = 0;
        prev = nullptr;
        current = this;

        YBI_ERROR(base != nullptr, "Virtual Memory Reservation Failed");
    }

    ~Arena()
    {
        if (current && current != this)
        {
            Arena *curr = current;
            while (curr != this)
            {
                Arena *toDelete = curr;
                curr = curr->prev;
                toDelete->prev = nullptr;
                toDelete->~Arena();
                AlignedFree(toDelete);
            }
            current = this;
            used = 0;
        }
        if (base)
        {
            Free(base, reserveSize, commits.data(), commits.size());
            base = nullptr;
        }
    }

    Arena(const Arena &) = delete;
    Arena &operator=(const Arena &) = delete;
    Arena(Arena &&) = delete;
    Arena &operator=(Arena &&) = delete;

    template <typename T>
    View<T> Push()
    {
        return PushArray(1);
    }

    template <typename T>
    View<T> PushArray(size_t count, size_t alignment = 0)
    {
        alignment = alignment == 0 ? current->alignment : alignment;
        size_t bytesNeeded = count * sizeof(T);

        size_t alignedUsed = AlignUp(current->used, alignment);
        size_t newUsed = alignedUsed + bytesNeeded;

        if (newUsed > current->reserveSize)
        {
            size_t newBlockSize = current->reserveSize;
            if (bytesNeeded > current->reserveSize)
            {
                newBlockSize = bytesNeeded;
            }

            void *newPtr = AlignedAlloc(sizeof(Arena), 8);
            Arena *newArena = new (newPtr) Arena(newBlockSize);
            newArena->prev = current;
            current = newArena;

            alignedUsed = 0;
            newUsed = bytesNeeded;
        }

        if (newUsed > current->committed)
        {
            size_t commitTarget =
                std::min(AlignUp(newUsed, current->commitSize), current->reserveSize);
            size_t commitRequestSize = commitTarget - current->committed;

            void *cookie = Commit(current->base + current->committed, commitRequestSize);
            current->commits.EmplaceBack(
                current->base + current->committed, commitRequestSize, cookie);
            current->committed = commitTarget;
        }

        T *dest = (T *)(current->base + alignedUsed);

        current->used = newUsed;
        return {dest, count};
    }
    void Clear()
    {
        Arena *curr = current;
        while (curr != this)
        {
            Arena *toDelete = curr;
            curr = curr->prev;
            toDelete->prev = nullptr;
            toDelete->~Arena();
            AlignedFree(toDelete);
        }
        current = this;
        used = 0;
    }

private:
    uint8_t *Reserve(size_t size)
    {
#if defined(_WIN32)
        return (uint8_t *)VirtualAlloc(nullptr, size, MEM_RESERVE, PAGE_NOACCESS);
#else
        return (uint8_t *)mmap(nullptr, size, PROT_NONE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
#endif
    }
    void *Commit(uint8_t *ptr, size_t size)
    {
#if defined(_WIN32)
        VirtualAlloc(ptr, size, MEM_COMMIT, PAGE_READWRITE);
#else
        mprotect(ptr, size, PROT_READ | PROT_WRITE);
#endif
        return nullptr;
    }
    void Free(uint8_t *base, size_t reserveSize, const ArenaCommit *commits, size_t numCommits)
    {
        (void)commits;
        (void)numCommits;
#if defined(_WIN32)
        VirtualFree(base, 0, MEM_RELEASE);
#else
        munmap(base, reserveSize);
#endif
    }
    size_t GetPageSize()
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

} // namespace ybi

#undef ARENA_RESERVE_SIZE
#undef ARENA_COMMIT_SIZE
