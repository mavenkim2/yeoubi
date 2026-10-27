#pragma once

#include "core/assert.h"
#ifdef YBI_ENABLE_CUDA

#define ARENA_RESERVE_SIZE (64ll * 1024ll * 1024ll)
#define ARENA_COMMIT_SIZE (64ll * 1024ll)

#include "core/bit.h"
#include "device/cuda_assert.h"
#include "device/device_memory.h"
#include <cstdint>
#include <cuda.h>
#include <cuda_runtime_api.h>

namespace ybi
{

class CUDAArena : public DeviceArena
{
private:
    struct ArenaCommit
    {
        uint8_t *ptr;
        size_t size;
        void *cookie;

        ArenaCommit(uint8_t *ptr, size_t size, void *cookie) : ptr(ptr), size(size), cookie(cookie)
        {
        }
    };

    CUDAArena *current;
    CUDAArena *prev;
    uint8_t *base;
    size_t used;
    size_t committed;
    size_t reserveSize;
    size_t alignment;
    size_t commitSize;
    Array<CUDAArena::ArenaCommit> commits;
    CUmemLocationType memType;

public:
    explicit CUDAArena(CUmemLocationType type = CU_MEM_LOCATION_TYPE_DEVICE,
                       size_t requestedReserveSize = ARENA_RESERVE_SIZE)
        : memType(type)
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

    ~CUDAArena()
    {
        if (current && current != this)
        {
            CUDAArena *curr = current;
            while (curr != this)
            {
                CUDAArena *toDelete = curr;
                curr = curr->prev;
                toDelete->prev = nullptr;
                toDelete->~CUDAArena();
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

    CUDAArena(const CUDAArena &) = delete;
    CUDAArena &operator=(const CUDAArena &) = delete;
    CUDAArena(CUDAArena &&) = delete;
    CUDAArena &operator=(CUDAArena &&) = delete;

    template <typename T>
    DeviceBuffer<T> Push()
    {
        return PushArray<T>(1);
    }

    template <typename T>
    DeviceBuffer<T> PushArray(size_t count, size_t alignment = 0)
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

            void *newPtr = AlignedAlloc(sizeof(CUDAArena), 8);
            CUDAArena *newArena = new (newPtr) CUDAArena(current->memType, newBlockSize);
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
        return DeviceBuffer<T>(dest, count);
    }
    void Clear()
    {
        CUDAArena *curr = current;
        while (curr != this)
        {
            CUDAArena *toDelete = curr;
            curr = curr->prev;
            toDelete->prev = nullptr;
            toDelete->~CUDAArena();
            AlignedFree(toDelete);
        }
        current = this;
        used = 0;
    }

protected:
    DeviceArena *GetCurrent() const override
    {
        return current;
    }

    size_t GetUsed() const override
    {
        return used;
    }

    void RestoreMark(Mark mark) override
    {
        YBI_ASSERT(this == mark.arena);
        CUDAArena *curr = current;

        while (curr != mark.current)
        {
            CUDAArena *toDelete = curr;
            curr = curr->prev;
            toDelete->prev = nullptr;
            toDelete->~CUDAArena();
            AlignedFree(toDelete);
        }
        current = (CUDAArena *)mark.current;
        current->used = mark.used;
    }

private:
    uint8_t *Reserve(size_t size)
    {
        CUdeviceptr dptr;
        CUDA_ASSERT(cuMemAddressReserve(&dptr, size, 0, 0, 0));
        return (uint8_t *)dptr;
    }

    void *Commit(uint8_t *ptr, size_t size)
    {
        CUmemGenericAllocationHandle handle;
        CUmemAllocationProp prop = {};
        prop.type = CU_MEM_ALLOCATION_TYPE_PINNED;
        prop.location.type = memType;
        prop.location.id = 0;

        CUDA_ASSERT(cuMemCreate(&handle, size, &prop, 0));
        CUDA_ASSERT(cuMemMap((CUdeviceptr)ptr, size, 0, handle, 0));
        CUmemAccessDesc access = {};
        access.location.type = CU_MEM_LOCATION_TYPE_DEVICE;
        access.location.id = 0;
        access.flags = CU_MEM_ACCESS_FLAGS_PROT_READWRITE;
        CUDA_ASSERT(cuMemSetAccess((CUdeviceptr)ptr, size, &access, 1));
        return (void *)(uintptr_t)handle;
    }

    void Free(uint8_t *base,
              size_t reserveSize,
              const CUDAArena::ArenaCommit *commits,
              size_t numCommits)
    {
        for (size_t i = 0; i < numCommits; i++)
        {
            CUDA_ASSERT(cuMemUnmap((CUdeviceptr)commits[i].ptr, commits[i].size));
            CUDA_ASSERT(cuMemRelease((CUmemGenericAllocationHandle)(uintptr_t)commits[i].cookie));
        }
        CUDA_ASSERT(cuMemAddressFree((CUdeviceptr)base, reserveSize));
    }

    size_t GetPageSize()
    {
        CUmemAllocationProp prop = {};
        prop.type = CU_MEM_ALLOCATION_TYPE_PINNED;
        prop.location.type = CU_MEM_LOCATION_TYPE_DEVICE;
        prop.location.id = 0;

        size_t granularity = 0;
        CUresult res =
            cuMemGetAllocationGranularity(&granularity, &prop, CU_MEM_ALLOC_GRANULARITY_MINIMUM);

        if (res != CUDA_SUCCESS)
        {
            return 0;
        }
        return granularity;
    }
};

} // namespace ybi

#undef ARENA_RESERVE_SIZE
#undef ARENA_COMMIT_SIZE

#endif
