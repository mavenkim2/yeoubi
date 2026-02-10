#pragma once

#include <cassert>
#include <cstddef>

YBI_NAMESPACE_BEGIN

template <typename T>
struct MemoryView
{
    T *ptr = nullptr;
    size_t count = 0;

    T &operator[](size_t index)
    {
        assert(index < count && "View access out of bounds!");
        return ptr[index];
    }
    T &operator[](size_t index) const
    {
        assert(index < count && "View access out of bounds!");
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
};

YBI_NAMESPACE_END
