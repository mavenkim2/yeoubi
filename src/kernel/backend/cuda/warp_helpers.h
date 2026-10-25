#pragma once

#include "core/vec3f.h"
#include "device/cuda_compat.h"

#if defined(__CUDACC__)

#include <type_traits>

namespace ybi
{

YBI_GPU uint32_t WarpInclusivePrefixSum(uint32_t value)
{
    const uint32_t threadIndex =
        threadIdx.x + blockDim.x * (threadIdx.y + blockDim.y * threadIdx.z);
    const uint32_t laneIndex = threadIndex & 31u;
#pragma unroll
    for (uint32_t offset = 1; offset < 32; offset <<= 1)
    {
        const uint32_t other = __shfl_up_sync(0xffffffffu, value, offset);
        if (laneIndex >= offset)
            value += other;
    }
    return value;
}

YBI_GPU uint32_t WarpExclusivePrefixSum(uint32_t value)
{
    return WarpInclusivePrefixSum(value) - value;
}

template <typename T>
YBI_GPU T WarpReduceSum(T value)
{
    static_assert(std::is_trivially_copyable<T>::value, "Warp reduction requires a trivially copyable type");
#pragma unroll
    for (uint32_t offset = 16; offset > 0; offset >>= 1)
    {
        uint32_t words[(sizeof(T) + sizeof(uint32_t) - 1) / sizeof(uint32_t)] = {};
        memcpy(words, &value, sizeof(T));
#pragma unroll
        for (uint32_t word = 0; word < sizeof(words) / sizeof(words[0]); ++word)
            words[word] = __shfl_down_sync(0xffffffffu, words[word], offset);
        T other = value;
        memcpy(&other, words, sizeof(T));
        value = value + other;
    }
    return value;
}

YBI_GPU float WarpReduceSum(float value)
{
#pragma unroll
    for (uint32_t offset = 16; offset > 0; offset >>= 1)
        value += __shfl_down_sync(0xffffffffu, value, offset);
    return value;
}

YBI_GPU float WarpReduceMin(float value)
{
#pragma unroll
    for (uint32_t offset = 16; offset > 0; offset >>= 1)
        value = fminf(value, __shfl_down_sync(0xffffffffu, value, offset));
    return value;
}

YBI_GPU float WarpReduceMax(float value)
{
#pragma unroll
    for (uint32_t offset = 16; offset > 0; offset >>= 1)
        value = fmaxf(value, __shfl_down_sync(0xffffffffu, value, offset));
    return value;
}

YBI_GPU Vec3f WarpReduceSum(Vec3f value)
{
    value.x = WarpReduceSum(value.x);
    value.y = WarpReduceSum(value.y);
    value.z = WarpReduceSum(value.z);
    return value;
}

YBI_GPU float WarpReduceSum(float value, uint32_t mask)
{
    const uint32_t threadIndex =
        threadIdx.x + blockDim.x * (threadIdx.y + blockDim.y * threadIdx.z);
    const uint32_t laneIndex = threadIndex & 31u;

#pragma unroll
    for (uint32_t offset = 16; offset > 0; offset >>= 1)
    {
        const float other = __shfl_down_sync(mask, value, offset);
        const uint32_t sourceLane = laneIndex + offset;
        if (sourceLane < 32 && (mask & (1u << sourceLane)))
            value += other;
    }
    return value;
}

YBI_GPU Vec3f WarpReduceSum(Vec3f value, uint32_t mask)
{
    value.x = WarpReduceSum(value.x, mask);
    value.y = WarpReduceSum(value.y, mask);
    value.z = WarpReduceSum(value.z, mask);
    return value;
}

YBI_GPU int64_t WarpReduceSum(int64_t value)
{
#pragma unroll
    for (uint32_t offset = 16; offset > 0; offset >>= 1)
        value += __shfl_down_sync(0xffffffffu, value, offset);
    return value;
}

} // namespace ybi

#endif
