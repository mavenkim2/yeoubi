#pragma once

#include "core/arena.h"

namespace ybi
{

struct ThreadContext
{
    Arena arena;
};

inline ThreadContext &GetThreadContext()
{
    thread_local ThreadContext context;
    return context;
}

struct ScratchArena
{
    Arena *arena;
    Arena::Scope scope;

    ScratchArena() : arena(&GetThreadContext().arena), scope(arena) {}
    ~ScratchArena() = default;

    ScratchArena(const ScratchArena &) = delete;
    ScratchArena &operator=(const ScratchArena &) = delete;
};

} // namespace ybi
