#pragma once

#ifdef WITH_CUDA

#include "cuda.h"
#include "cuda_runtime_api.h"
#include "device/cuda_assert.h"
#include "device/device_memory_view.h"
#include "util/memory_arena.h"

YBI_NAMESPACE_BEGIN

struct CUDAMemoryProvider
{
    template <typename T>
    using View = DeviceMemoryView<T>;

    static uint8_t *Reserve(size_t size)
    {
        CUdeviceptr dptr;
        CUDA_ASSERT(cuMemAddressReserve(&dptr, size, 0, 0, 0));
        return (uint8_t *)dptr;
    }

    static void Commit(uint8_t *ptr, size_t size)
    {
        CUmemGenericAllocationHandle handle;
        CUmemAllocationProp prop = {};
        prop.type = CU_MEM_ALLOCATION_TYPE_PINNED;
        prop.location.type = CU_MEM_LOCATION_TYPE_DEVICE;
        prop.location.id = 0;

        CUDA_ASSERT(cuMemCreate(&handle, size, &prop, 0));
        CUDA_ASSERT(cuMemMap((CUdeviceptr)ptr, size, 0, handle, 0));
        CUmemAccessDesc access = {};
        access.location.type = CU_MEM_LOCATION_TYPE_DEVICE;
        access.location.id = 0;
        access.flags = CU_MEM_ACCESS_FLAGS_PROT_READWRITE;
        CUDA_ASSERT(cuMemSetAccess((CUdeviceptr)ptr, size, &access, 1));
    }

    static void Free(uint8_t *ptr, size_t size)
    {
        CUDA_ASSERT(cuMemUnmap((CUdeviceptr)ptr, size));
        CUDA_ASSERT(cuMemAddressFree((CUdeviceptr)ptr, size));
    }

    static size_t GetPageSize()
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

using CUDAMemoryArena = MemoryArena<CUDAMemoryProvider>;

YBI_NAMESPACE_END

#endif
