#pragma once

#include "core/arena.h"
#include "core/assert.h"
#include "core/bit.h"
#include "core/math.h"
#include <cstdint>
#include <cstring>
#include <mutex>
#include <thread>

namespace ybi
{

template <typename T>
struct RingBuffer
{
    View<T> entries;
    uint64_t readOffset;
    uint64_t writeOffset;

    uint32_t maxEntries;

    RingBuffer(Arena *arena, uint32_t maxEntries);
    bool Write(T *vals, uint32_t num);
    void WriteWithOverwrite(T *vals, uint32_t num);
    void SynchronizedWrite(std::mutex *mutex, T *vals, uint32_t num);
    T *Read(Arena *arena, uint32_t &num);
    T *SynchronizedRead(std::mutex *mutex, Arena *arena, uint32_t &num);
};


template <typename T>
RingBuffer<T>::RingBuffer(Arena *arena, uint32_t maxEntries) : maxEntries(maxEntries)
{
    YBI_ASSERT(IsPowerOfTwo(maxEntries));
    entries = arena->PushArray<T>(maxEntries);
    readOffset  = 0;
    writeOffset = 0;
}

template <typename T>
bool RingBuffer<T>::Write(T *vals, uint32_t num)
{
    bool result  = true;
    uint32_t capacity = static_cast<uint32_t>(entries.count);
    if (writeOffset + num >= readOffset + capacity)
    {
        result = false;
    }
    else
    {
        uint32_t writeIndex = static_cast<uint32_t>(writeOffset & (capacity - 1));
        uint32_t numToEnd   = capacity - writeIndex;
        memcpy(entries.data() + writeIndex, vals, sizeof(T) * Min(num, numToEnd));
        if (num > numToEnd)
        {
            memcpy(entries.data(), vals + numToEnd, sizeof(T) * (num - numToEnd));
        }
        writeOffset += num;
    }

    return result;
}

template <typename T>
void RingBuffer<T>::WriteWithOverwrite(T *vals, uint32_t num)
{
    uint32_t capacity = static_cast<uint32_t>(entries.count);
    if (writeOffset + num >= readOffset + capacity)
    {
        readOffset = writeOffset + num - capacity;
    }
    uint32_t writeIndex = static_cast<uint32_t>(writeOffset & (capacity - 1));
    uint32_t numToEnd   = capacity - writeIndex;
    memcpy(entries.data() + writeIndex, vals, sizeof(T) * Min(num, numToEnd));
    if (num > numToEnd)
    {
        memcpy(entries.data(), vals + numToEnd, sizeof(T) * (num - numToEnd));
    }
    writeOffset += num;
}

template <typename T>
void RingBuffer<T>::SynchronizedWrite(std::mutex *mutex, T *vals, uint32_t num)
{
    for (;;)
    {
        {
            std::lock_guard<std::mutex> lock(*mutex);
            if (Write(vals, num))
            {
                break;
            }
        }

        std::this_thread::yield();
    }
}

template <typename T>
T *RingBuffer<T>::Read(Arena *arena, uint32_t &num)
{
    YBI_ASSERT(readOffset <= writeOffset);
    uint32_t numToRead = static_cast<uint32_t>(writeOffset - readOffset);
    numToRead     = num == 0 ? numToRead : Min(numToRead, num);
    uint32_t capacity  = static_cast<uint32_t>(entries.count);
    T *vals       = 0;
    if (numToRead)
    {
        vals          = arena->PushArray<T>(numToRead).data();
        uint32_t readIndex = static_cast<uint32_t>(readOffset & (capacity - 1));
        uint32_t numToEnd  = capacity - readIndex;
        memcpy(vals, entries.data() + readIndex, sizeof(T) * Min(numToRead, numToEnd));
        if (numToRead > numToEnd)
        {
            memcpy(vals + numToEnd, entries.data(), sizeof(T) * (numToRead - numToEnd));
        }
        readOffset += numToRead;
    }
    num = numToRead;
    return vals;
}

template <typename T>
T *RingBuffer<T>::SynchronizedRead(std::mutex *mutex, Arena *arena, uint32_t &num)
{
    std::lock_guard<std::mutex> lock(*mutex);
    T *result = Read(arena, num);
    return result;
}


} // namespace ybi
