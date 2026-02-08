#pragma once

#include "aligned_malloc.h"
#include <assert.h>
#include <type_traits>

YBI_NAMESPACE_BEGIN

template <typename C, typename = void>
struct has_data_and_size : std::false_type
{
};

template <typename C>
struct has_data_and_size<
    C,
    std::void_t<decltype(std::declval<C>().data()), decltype(std::declval<C>().size())>>
    : std::true_type
{
};

template <typename T, size_t alignment = 16>
class Array
{
    using value_type = T;
    static_assert((alignment & (alignment - 1)) == 0, "Alignment must be a power of 2");

public:
    Array() : m_data(nullptr), m_size(0), m_capacity(0) {}
    explicit Array(size_t requestedSize) : m_data(nullptr), m_size(0), m_capacity(0)
    {
        if (requestedSize == 0)
        {
            m_data = nullptr;
            m_size = 0;
            m_capacity = 0;
        }
        else
        {
            m_data = AllocateAligned(requestedSize);
            m_size = requestedSize;
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

    Array(Array &&other) : m_data(other.m_data), m_size(other.m_size), m_capacity(other.m_capacity)
    {
        other.m_data = nullptr;
        other.m_size = 0;
        other.m_capacity = 0;
    }

    template <typename OtherContainer,
              typename = std::enable_if_t<has_data_and_size<OtherContainer>::value>>
    Array(const OtherContainer &other) : m_data(nullptr), m_size(0), m_capacity(0)
    {
        CopyFrom(other);
    }

    // Templated assignment operator
    template <typename OtherContainer,
              typename = std::enable_if_t<has_data_and_size<OtherContainer>::value>>
    Array &operator=(const OtherContainer &other)
    {
        CopyFrom(other);
        return *this;
    }

    void Resize(size_t newSize)
    {
        if (newSize > m_capacity)
        {
            size_t newCapacity = newSize;
            T *newBuffer = AllocateAligned(newCapacity);
            assert(newBuffer);

            if (m_data)
            {
                memcpy(newBuffer, m_data, m_size * sizeof(T));
                util::AlignedFree(m_data);
            }

            m_data = newBuffer;
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

            m_data = newBuffer;
            m_capacity = newCapacity;
        }
    }

    T &operator[](size_t index)
    {
        assert(index < m_size);
        return m_data[index];
    }
    const T &operator[](size_t index) const
    {
        assert(index < m_size);
        return m_data[index];
    }
    T *data()
    {
        return m_data;
    }
    const T *data() const
    {
        return m_data;
    }
    size_t size() const
    {
        return m_size;
    }
    size_t capacity() const
    {
        return m_capacity;
    }
    T *begin()
    {
        return m_data;
    }
    const T *begin() const
    {
        return m_data;
    }
    T *end()
    {
        return m_data + m_size;
    }
    const T *end() const
    {
        return m_data + m_size;
    }

private:
    // TODO: investigate performance
    T *AllocateAligned(size_t n)
    {
        assert(n != 0);
        void *ptr = nullptr;
        ptr = util::AlignedAlloc(sizeof(T) * n, alignment);
        assert(ptr);
        return static_cast<T *>(ptr);
    }

    template <typename OtherContainer>
    void CopyFrom(const OtherContainer &other)
    {
        static_assert(sizeof(typename OtherContainer::value_type) == sizeof(T),
                      "Container element sizes must match");

        size_t newCount = other.size();
        Resize(newCount);
        if (newCount > 0)
        {
            memcpy(m_data, other.data(), newCount * sizeof(T));
        }
    }

    T *m_data;
    size_t m_size;
    size_t m_capacity;
};

YBI_NAMESPACE_END
