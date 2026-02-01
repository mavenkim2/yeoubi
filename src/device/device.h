#pragma once
#include <cstdint>
#include <string>

YBI_NAMESPACE_BEGIN

struct BVH;

typedef uint32_t ModuleHandle;
typedef uint32_t KernelHandle;

struct GPUArena
{
    virtual void *Alloc(size_t size, uintptr_t alignment) = 0;
    virtual void Clear()                                  = 0;

    template <typename T>
    T *Alloc(uint32_t count, uint32_t alignment = 0)
    {
        alignment = alignment == 0 ? sizeof(T) : alignment;
        return (T *)Alloc(sizeof(T) * count, alignment);
    }
};

struct Device
{
    virtual ModuleHandle RegisterModule(std::string module)                       = 0;
    virtual KernelHandle RegisterKernels(const char *kernel, ModuleHandle module) = 0;
    template <typename... Args>
    void ExecuteKernel(KernelHandle handle, uint32_t numBlocks, uint32_t blockSize,
                       Args... args)
    {
        void *paramArray[] = {(void *)&args...};
        ExecuteKernelInternal(handle, numBlocks, blockSize, paramArray, sizeof...(args));
    }

    virtual void *Alloc(uint32_t size, uintptr_t alignment) = 0;
    template <typename T>
    T *Alloc(uint32_t count, uint32_t alignment = 0)
    {
        alignment = alignment == 0 ? sizeof(T) : alignment;
        return (T *)Alloc(sizeof(T) * count, alignment);
    }
    virtual void MemZero(void *ptr, uint64_t size)         = 0;
    virtual void MemSet(void *ptr, char ch, uint64_t size) = 0;
    virtual void BuildBVH(BVH *bvh)                        = 0;
    virtual GPUArena *CreateArena(size_t maxSize)          = 0;

protected:
    virtual void ExecuteKernelInternal(KernelHandle handle, uint32_t numBlocks,
                                       uint32_t blockSize, void **params,
                                       uint32_t paramCount) = 0;
};

// struct Module
// {
//     string name;
//     string data;
// };
//
// struct Kernel
// {
//     std::string name;
//     ModuleHandle moduleHandle;
// };

YBI_NAMESPACE_END
