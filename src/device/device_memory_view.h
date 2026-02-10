#pragma once

#include "util/assert.h"
#include <cstddef>

YBI_NAMESPACE_BEGIN

template <typename T>
struct DeviceMemoryView
{
    using value_type = T;

    T *ptr = nullptr;
    size_t count = 0;

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
    DeviceMemoryView<T> operator+(size_t offset) const
    {
        ASSERT(offset <= count);

        DeviceMemoryView<T> result;
        result.ptr = ptr + offset;
        result.count = count - offset;
        return result;
    }
    DeviceMemoryView<T> &operator+=(size_t offset)
    {
        ASSERT(offset <= count);
        ptr += offset;
        count -= offset;
        return *this;
    }
};

YBI_NAMESPACE_END
