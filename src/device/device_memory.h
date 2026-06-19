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

struct DeviceBuffer
{
    DevicePtr<std::byte> ptr;
    size_t sizeInBytes = 0;

    constexpr DeviceBuffer() = default;
    constexpr DeviceBuffer(std::nullptr_t) : ptr(nullptr), sizeInBytes(0) {}
    DeviceBuffer(void *ptr, size_t sizeInBytes)
        : ptr(DevicePtr<std::byte>(static_cast<std::byte *>(ptr))), sizeInBytes(sizeInBytes)
    {
    }

    DeviceBuffer(uintptr_t address, size_t sizeInBytes)
        : ptr(DevicePtr<std::byte>(address)), sizeInBytes(sizeInBytes)
    {
    }

    std::byte *Data() const
    {
        return ptr.Get();
    }

    uintptr_t Address() const
    {
        return ptr.Address();
    }

    size_t SizeInBytes() const
    {
        return sizeInBytes;
    }

    bool Empty() const
    {
        return !ptr || sizeInBytes == 0;
    }

    explicit operator bool() const
    {
        return !Empty();
    }

    template <typename T>
    DevicePtr<T> As() const
    {
        static_assert(!std::is_void_v<T>, "Use Address() for untyped device addresses");
        YBI_ASSERT(Address() % alignof(T) == 0);
        YBI_ASSERT(sizeInBytes % sizeof(T) == 0);
        return ptr.As<T>();
    }

    template <typename T>
    size_t Count() const
    {
        static_assert(!std::is_void_v<T>, "void has no element count");
        YBI_ASSERT(sizeInBytes % sizeof(T) == 0);
        return sizeInBytes / sizeof(T);
    }
};

} // namespace ybi
