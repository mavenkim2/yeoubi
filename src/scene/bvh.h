#pragma once

#include <cstdint>
namespace ybi
{

struct BVHHandle
{
    uint64_t handle;
    BVHHandle()
    {
        handle = ~0ull;
    }
    BVHHandle(uint64_t handle) : handle(handle) {}
    bool IsValid()
    {
        return handle != ~0ull;
    }
};

struct BVH
{
    BVHHandle handle;

    BVH(BVHHandle handle) : handle(handle) {}
};

} // namespace ybi
