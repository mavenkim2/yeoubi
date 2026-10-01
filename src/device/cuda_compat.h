#pragma once

#if defined(__CUDACC_RTC__)
#include <cuda/std/cmath>
#include <cuda/std/cstddef>
#include <cuda/std/cstdint>
#include <cuda/std/cstring>

using cuda::std::memcpy;
using cuda::std::memmove;
using cuda::std::memset;
using cuda::std::ptrdiff_t;
using cuda::std::size_t;
using cuda::std::int8_t;
using cuda::std::int16_t;
using cuda::std::int32_t;
using cuda::std::int64_t;
using cuda::std::intptr_t;
using cuda::std::uint8_t;
using cuda::std::uint16_t;
using cuda::std::uint32_t;
using cuda::std::uint64_t;
using cuda::std::uintptr_t;
#else
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>

using std::memcpy;
using std::memmove;
using std::memset;
using std::ptrdiff_t;
using std::size_t;
using std::int8_t;
using std::int16_t;
using std::int32_t;
using std::int64_t;
using std::intptr_t;
using std::uint8_t;
using std::uint16_t;
using std::uint32_t;
using std::uint64_t;
using std::uintptr_t;
#endif

#ifndef YBI_DEVICE

#if defined(__CUDACC__)
#define YBI_GPU __device__ __forceinline__
#define YBI_CPU_GPU __host__ __device__ __forceinline__
#define YBI_TABLE_GPU __device__
#else
#define YBI_GPU inline
#define YBI_CPU_GPU inline
#define YBI_TABLE_GPU
#endif

#endif
