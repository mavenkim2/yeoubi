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
