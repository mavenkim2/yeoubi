#pragma once

#include "core/arena.h"
#include "core/assert.h"
#include "core/math.h"
#include <mutex>
#include <vector>

namespace ybi
{

struct ConcurrentHashMap
{
    struct Node
    {
        size_t index;
        Node *next;
    };

    View<Node *> nodes;
    std::vector<std::mutex> mutexes;

    ConcurrentHashMap(Arena *arena, uint32_t numNodes, uint32_t numShards = 0);

    template <typename CompareFn>
    size_t GetOrAdd(Arena *arena, size_t index, uint64_t hash, CompareFn compareFn);

    template <typename CompareFn, typename AddFn>
    size_t GetOrAdd(Arena *arena, size_t index, uint64_t hash, CompareFn compareFn, AddFn addFn);
};

inline ConcurrentHashMap::ConcurrentHashMap(Arena *arena, uint32_t numNodes, uint32_t inputNumShards)
    : mutexes(inputNumShards == 0 ? Max(1u, numNodes / 8) : inputNumShards)
{
    nodes = arena->PushArrayZero<Node *>(numNodes);
}

template <typename CompareFn>
size_t ConcurrentHashMap::GetOrAdd(Arena *arena, size_t index, uint64_t hash, CompareFn compareFn)
{
    return GetOrAdd(arena, index, hash, compareFn, [](size_t) {});
}

template <typename CompareFn, typename AddFn>
size_t ConcurrentHashMap::GetOrAdd(Arena *arena,
                                   size_t index,
                                   uint64_t hash,
                                   CompareFn compareFn,
                                   AddFn addFn)
{
    std::mutex &mutex = mutexes[hash % mutexes.size()];

    std::lock_guard<std::mutex> lock(mutex);

    Node **ptr = &nodes[hash % nodes.size()];
    Node *node = *ptr;
    while (node)
    {
        if (compareFn(node->index))
        {
            return node->index;
        }
        ptr = &node->next;
        node = node->next;
    }

    YBI_ASSERT(node == nullptr);
    addFn(index);

    Node *nextNode = arena->PushArrayZero<Node>(1).data();
    nextNode->index = index;
    *ptr = nextNode;

    return index;
}

} // namespace ybi
