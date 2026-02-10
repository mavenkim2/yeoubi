#pragma once

#include "util/aligned_malloc.h"
#include "util/assert.h"
#include "util/base.h"
#include <algorithm>
#include <cstdint>
#include <cstring>

#define ARENA_RESERVE_SIZE (64ll * 1024ll * 1024ll)
#define ARENA_COMMIT_SIZE (64ll * 1024ll)

YBI_NAMESPACE_BEGIN

template <typename Provider>
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
    template <typename T>
    using ViewType = typename Provider::template View<T>;

    explicit MemoryArena(size_t requestedReserveSize = ARENA_RESERVE_SIZE)
    {
        size_t pageSize = Provider::GetPageSize();
        reserveSize = util::AlignUp(requestedReserveSize, pageSize);
        commitSize = util::AlignUp(ARENA_COMMIT_SIZE, pageSize);

        ERROR(reserveSize % commitSize == 0, "Reserve size should be multiple of commit size.\n");
        alignment = 16;

        base = Provider::Reserve(reserveSize);
        used = 0;
        committed = 0;
        prev = nullptr;
        current = this;

        ERROR(base != nullptr, "Virtual Memory Reservation Failed");
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
            Provider::Free(base, reserveSize);
        }
    }

    MemoryArena(const MemoryArena &) = delete;
    MemoryArena &operator=(const MemoryArena &) = delete;
    MemoryArena(MemoryArena &&) = delete;
    MemoryArena &operator=(MemoryArena &&) = delete;

    template <typename T>
    ViewType<T> Push()
    {
        return PushArray(1);
    }

    ViewType<uint8_t> PushBytes(size_t numBytes)
    {
        return PushArray<uint8_t>(numBytes);
    }

    template <typename T>
    ViewType<T> PushArray(size_t count)
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

            Provider::Commit(current->base + current->committed, commitRequestSize);
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
