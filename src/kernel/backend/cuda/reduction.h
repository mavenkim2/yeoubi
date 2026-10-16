#pragma once

#include "kernel/backend/cuda/warp_reduce.h"

#if defined(__CUDACC__)

namespace ybi
{

// T{} must be the identity for operator+. Only thread zero receives the block sum.
template <typename T, uint32_t BlockSize = 256>
YBI_GPU T BlockReduceSum(T value)
{
    static_assert(BlockSize >= 32 && BlockSize <= 1024 && BlockSize % 32 == 0,
                  "Block reduction requires complete warps");
    static_assert(std::is_trivially_copyable<T>::value, "Block reduction requires a trivially copyable type");
    static_assert(std::is_trivially_default_constructible<T>::value,
                  "Shared reduction storage requires a trivial default constructor");
    YBI_ASSERT(blockDim.x == BlockSize && blockDim.y == 1 && blockDim.z == 1);
    __shared__ T warpValues[BlockSize / 32];
    const uint32_t lane = threadIdx.x & 31u;
    const uint32_t warp = threadIdx.x >> 5;
    value = WarpReduceSum(value);
    if (lane == 0)
        warpValues[warp] = value;
    __syncthreads();

    T result = {};
    if (warp == 0)
    {
        if (lane < BlockSize / 32)
            result = warpValues[lane];
        result = WarpReduceSum(result);
    }
    __syncthreads();
    return result;
}

// Each pass writes CeilDiv(count, BlockSize) values to a separate output buffer.
template <typename T, uint32_t BlockSize = 256>
__global__ void ReduceKernel(const T *input, uint32_t count, T *output)
{
    const uint32_t outputCount = CeilDiv(count, BlockSize);
    for (uint32_t block = blockIdx.x; block < outputCount; block += gridDim.x)
    {
        const uint32_t index = block * BlockSize + threadIdx.x;
        T value = {};
        if (index < count)
            value = input[index];
        const T sum = BlockReduceSum<T, BlockSize>(value);
        if (threadIdx.x == 0)
            output[block] = sum;
    }
}

} // namespace ybi

#endif
