#pragma once

#include "device/cuda_compat.h"

#if defined(__CUDACC__)
#include <cuda/atomic>
#endif

#if defined(_MSC_VER) && !defined(__CUDACC_RTC__)
#include <intrin.h>
#endif

namespace ybi
{

enum class AtomicOrder
{
    Relaxed,
    Acquire,
    Release,
    AcqRel,
    SeqCst,
};

#if defined(__CUDA_ARCH__)
YBI_GPU cuda::std::memory_order NativeAtomicOrder(AtomicOrder order)
{
    switch (order)
    {
        case AtomicOrder::Relaxed:
            return cuda::std::memory_order_relaxed;
        case AtomicOrder::Acquire:
            return cuda::std::memory_order_acquire;
        case AtomicOrder::Release:
            return cuda::std::memory_order_release;
        case AtomicOrder::AcqRel:
            return cuda::std::memory_order_acq_rel;
        case AtomicOrder::SeqCst:
            return cuda::std::memory_order_seq_cst;
    }
    return cuda::std::memory_order_seq_cst;
}
#elif !defined(_MSC_VER)
YBI_CPU_GPU int NativeAtomicOrder(AtomicOrder order)
{
    switch (order)
    {
        case AtomicOrder::Relaxed:
            return __ATOMIC_RELAXED;
        case AtomicOrder::Acquire:
            return __ATOMIC_ACQUIRE;
        case AtomicOrder::Release:
            return __ATOMIC_RELEASE;
        case AtomicOrder::AcqRel:
            return __ATOMIC_ACQ_REL;
        case AtomicOrder::SeqCst:
            return __ATOMIC_SEQ_CST;
    }
    return __ATOMIC_SEQ_CST;
}
#endif

// CUDA operations synchronize threads on the same device. MSVC's Interlocked
// operations provide a full barrier, which can be stronger than the requested order.
// All integer operations return the value before modification.
YBI_CPU_GPU uint32_t AtomicOr(uint32_t *destination,
                              uint32_t value,
                              AtomicOrder order = AtomicOrder::Relaxed)
{
#if defined(__CUDA_ARCH__)
    return cuda::atomic_ref<uint32_t, cuda::thread_scope_device>(*destination)
        .fetch_or(value, NativeAtomicOrder(order));
#elif defined(_MSC_VER)
    (void)order;
    return uint32_t(
        _InterlockedOr(reinterpret_cast<volatile long *>(destination), static_cast<long>(value)));
#else
    return __atomic_fetch_or(destination, value, NativeAtomicOrder(order));
#endif
}

YBI_CPU_GPU uint32_t AtomicAdd(uint32_t *destination,
                               uint32_t value,
                               AtomicOrder order = AtomicOrder::Relaxed)
{
#if defined(__CUDA_ARCH__)
    return cuda::atomic_ref<uint32_t, cuda::thread_scope_device>(*destination)
        .fetch_add(value, NativeAtomicOrder(order));
#elif defined(_MSC_VER)
    (void)order;
    return uint32_t(_InterlockedExchangeAdd(reinterpret_cast<volatile long *>(destination),
                                            static_cast<long>(value)));
#else
    return __atomic_fetch_add(destination, value, NativeAtomicOrder(order));
#endif
}

YBI_CPU_GPU uint32_t AtomicSub(uint32_t *destination,
                               uint32_t value,
                               AtomicOrder order = AtomicOrder::Relaxed)
{
    return AtomicAdd(destination, 0u - value, order);
}

YBI_CPU_GPU uint32_t AtomicCompareExchange(uint32_t *destination,
                                           uint32_t expected,
                                           uint32_t desired,
                                           AtomicOrder order = AtomicOrder::Relaxed)
{
#if defined(__CUDA_ARCH__)
    cuda::atomic_ref<uint32_t, cuda::thread_scope_device>(*destination)
        .compare_exchange_strong(expected, desired, NativeAtomicOrder(order));
    return expected;
#elif defined(_MSC_VER)
    (void)order;
    return uint32_t(_InterlockedCompareExchange(reinterpret_cast<volatile long *>(destination),
                                                static_cast<long>(desired),
                                                static_cast<long>(expected)));
#else
    const AtomicOrder failureOrder = order == AtomicOrder::Release  ? AtomicOrder::Relaxed
                                     : order == AtomicOrder::AcqRel ? AtomicOrder::Acquire
                                                                    : order;
    __atomic_compare_exchange_n(destination,
                                &expected,
                                desired,
                                false,
                                NativeAtomicOrder(order),
                                NativeAtomicOrder(failureOrder));
    return expected;
#endif
}

YBI_CPU_GPU float
AtomicAdd(float *destination, float value, AtomicOrder order = AtomicOrder::Relaxed)
{
#if defined(__CUDA_ARCH__)
    if (order == AtomicOrder::Relaxed)
        return atomicAdd(destination, value);
    return cuda::atomic_ref<float, cuda::thread_scope_device>(*destination)
        .fetch_add(value, NativeAtomicOrder(order));
#elif defined(_MSC_VER)
    (void)order;
    auto *bits = reinterpret_cast<volatile long *>(destination);
    long expected = _InterlockedCompareExchange(bits, 0, 0);
    for (;;)
    {
        float previous;
        memcpy(&previous, &expected, sizeof(previous));
        const float next = previous + value;
        long desired;
        memcpy(&desired, &next, sizeof(desired));
        const long observed = _InterlockedCompareExchange(bits, desired, expected);
        // Compare the bit patterns so a stored NaN cannot make the loop spin forever.
        if (observed == expected)
            return previous;
        expected = observed;
    }
#else
    float expected;
    __atomic_load(destination, &expected, __ATOMIC_RELAXED);
    for (;;)
    {
        float desired = expected + value;
        if (__atomic_compare_exchange(destination,
                                      &expected,
                                      &desired,
                                      false,
                                      NativeAtomicOrder(order),
                                      __ATOMIC_RELAXED))
            return expected;
    }
#endif
}

#if defined(__CUDACC__)

YBI_GPU void AtomicMin(float *destination, float value)
{
    const int bits = __float_as_int(value);
    if (bits >= 0)
        atomicMin(reinterpret_cast<int *>(destination), bits);
    else
        atomicMax(reinterpret_cast<unsigned int *>(destination), __float_as_uint(value));
}

YBI_GPU void AtomicMax(float *destination, float value)
{
    const int bits = __float_as_int(value);
    if (bits >= 0)
        atomicMax(reinterpret_cast<int *>(destination), bits);
    else
        atomicMin(reinterpret_cast<unsigned int *>(destination), __float_as_uint(value));
}

#endif

} // namespace ybi
