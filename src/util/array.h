#pragma once

#include "aligned_malloc.h"
#include <assert.h>

template <typename T, size_t alignment = 16>
class Array
{
    static_assert((alignment & (alignment - 1)) == 0, "Alignment must be a power of 2");

public:
    Array() : m_data(nullptr), m_size(0), m_capacity(0) {}
    explicit Array(size_t requestedSize) : m_data(nullptr), m_size(0), m_capacity(0)
    {
        if (requestedSize == 0)
        {
            m_data     = nullptr;
            m_size     = 0;
            m_capacity = 0;
        }
        else
        {
            m_data     = AllocateAligned(requestedSize);
            m_size     = requestedSize;
            m_capacity = requestedSize;
        }
    }

    ~Array()
    {
        if (m_data)
        {
            util::AlignedFree(m_data);
        }
    }

    Array(const Array &)            = delete;
    Array &operator=(const Array &) = delete;

    Array(Array &&other)
        : m_data(other.m_data), m_size(other.m_size), m_capacity(other.m_capacity)
    {
        other.m_data     = nullptr;
        other.m_size     = 0;
        other.m_capacity = 0;
    }

    void Resize(size_t newSize)
    {
        if (newSize > m_capacity)
        {
            size_t newCapacity = newSize;
            T *newBuffer       = AllocateAligned(newCapacity);
            assert(newBuffer);

            if (m_data)
            {
                memcpy(newBuffer, m_data, m_size * sizeof(T));
                util::AlignedFree(m_data);
            }

            m_data     = newBuffer;
            m_capacity = newCapacity;
        }

        m_size = newSize;
    }

    void Reserve(size_t newCapacity)
    {
        if (newCapacity > m_capacity)
        {
            T *newBuffer = AllocateAligned(newCapacity);
            assert(newBuffer);

            if (m_data)
            {
                memcpy(newBuffer, m_data, m_size * sizeof(T));
                util::AlignedFree(m_data);
            }

            m_data     = newBuffer;
            m_capacity = newCapacity;
        }
    }

    T &operator[](size_t index) { return m_data[index]; }
    const T &operator[](size_t index) const { return m_data[index]; }
    T *data() { return m_data; }
    size_t size() const { return m_size; }
    size_t capacity() const { return m_capacity; }

private:
    // TODO: investigate performance
    T *AllocateAligned(size_t n)
    {
        assert(n != 0);
        void *ptr = nullptr;
        ptr       = util::AlignedAlloc(n, alignment);
        assert(ptr);
        return static_cast<T *>(ptr);
    }

    T *m_data;
    size_t m_size;
    size_t m_capacity;
};
