#pragma once

#include "core/assert.h"
#include <cstddef>
#include <cstdint>
#include <type_traits>

namespace ybi
{

template <typename T>
struct DevicePtr
{
    T *ptr = nullptr;

    constexpr DevicePtr() = default;
    constexpr DevicePtr(std::nullptr_t) : ptr(nullptr) {}
    explicit constexpr DevicePtr(T *ptr) : ptr(ptr) {}
    explicit DevicePtr(uintptr_t address) : ptr(reinterpret_cast<T *>(address)) {}

    T *Get() const
    {
        return ptr;
    }

    uintptr_t Address() const
    {
        return reinterpret_cast<uintptr_t>(ptr);
    }

    explicit operator bool() const
    {
        return ptr != nullptr;
    }

    template <typename U>
    DevicePtr<U> As() const
    {
        return DevicePtr<U>(reinterpret_cast<U *>(ptr));
    }
};

template <typename T>
struct DeviceBuffer
{
    static_assert(!std::is_void_v<T>, "Use std::byte for raw device storage");

    DevicePtr<T> ptr;
    size_t count = 0;

    constexpr DeviceBuffer() = default;
    constexpr DeviceBuffer(std::nullptr_t) : ptr(nullptr), count(0) {}
    DeviceBuffer(T *ptr, size_t count) : ptr(DevicePtr<T>(ptr)), count(count)
    {
    }

    DeviceBuffer(uintptr_t address, size_t count) : ptr(DevicePtr<T>(address)), count(count)
    {
    }

    T *Data()
    {
        return ptr.Get();
    }

    T *Data() const
    {
        return ptr.Get();
    }

    void *KernelArg()
    {
        return &ptr.ptr;
    }

    uintptr_t Address() const
    {
        return ptr.Address();
    }

    size_t SizeInBytes() const
    {
        return count * sizeof(T);
    }

    bool Empty() const
    {
        return !ptr || count == 0;
    }

    explicit operator bool() const
    {
        return !Empty();
    }

    template <typename U>
    DeviceBuffer<U> As() const
    {
        static_assert(!std::is_void_v<U>, "Use std::byte for raw device storage");
        YBI_ASSERT(Address() % alignof(U) == 0);
        YBI_ASSERT(SizeInBytes() % sizeof(U) == 0);
        return DeviceBuffer<U>(ptr.template As<U>().Get(), SizeInBytes() / sizeof(U));
    }

    size_t Count() const
    {
        return count;
    }
};

} // namespace ybi
