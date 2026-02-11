#pragma once

#include "util/assert.h"
#include <cstddef>

YBI_NAMESPACE_BEGIN

template <typename T>
struct MemoryView
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

    MemoryView<T> operator+(size_t offset) const
    {
        YBI_ASSERT(offset <= count);

        MemoryView<T> result;
        result.ptr = ptr + offset;
        result.count = count - offset;
        return result;
    }
    MemoryView<T> &operator+=(size_t offset)
    {
        YBI_ASSERT(offset <= count);
        ptr += offset;
        count -= offset;
        return *this;
    }
};

YBI_NAMESPACE_END
