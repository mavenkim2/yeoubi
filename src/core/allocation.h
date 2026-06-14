#pragma once

#include <cstddef>
#include <new>
#include <type_traits>
#include <utility>

namespace ybi
{

void *AlignedAlloc(size_t size, size_t alignment);
void AlignedFree(void *ptr);

template <typename T>
T *AllocateArray(size_t count, size_t alignment = alignof(T))
{
    if (count == 0)
    {
        return nullptr;
    }

    const size_t size = sizeof(T) * count;
    void *memory = AlignedAlloc(size, alignment);
    if (!memory)
    {
        return nullptr;
    }

    T *items = static_cast<T *>(memory);
    if constexpr (!std::is_trivially_default_constructible_v<T>)
    {
        size_t constructed = 0;
        for (; constructed < count; ++constructed)
        {
            new (&items[constructed]) T();
        }
    }

    return items;
}

template <typename T>
void DestroyArray(T *items, size_t count)
{
    if (!items)
    {
        return;
    }

    if constexpr (!std::is_trivially_destructible_v<T>)
    {
        for (size_t i = 0; i < count; ++i)
        {
            items[i].~T();
        }
    }

    AlignedFree(items);
}

} // namespace ybi
