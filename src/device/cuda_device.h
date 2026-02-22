#pragma once

#if defined(WITH_CUDA) && defined(WITH_OPTIX)

#include "device/cuda_assert.h"
#include "device/cuda_device_memory_arena.h"
#include "scene/clusterizer.h"
#include "scene/micropolygon_mesh.h"
#include "util/aligned_malloc.h"
#include "util/array.h"
#include "util/assert.h"
#include "util/base.h"
#include "util/host_memory_arena.h"
#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#include <memory>
#include <optix_stubs.h>
#include <optix_types.h>
#include <string>

YBI_NAMESPACE_BEGIN

struct Scene;
struct ScenePool;
struct Mesh;
struct Curves;

struct OptixPrimaryPipelineState
{
    OptixPipeline pipeline = nullptr;
    OptixShaderBindingTable sbt = {};
    CUdeviceptr raygenRecordBuffer = 0;
    CUdeviceptr missRecordBuffer = 0;
    CUdeviceptr hitgroupRecordBuffer = 0;
    OptixModule module = nullptr;
    OptixModule curveModule = nullptr;
    OptixProgramGroup raygenGroup = nullptr;
    OptixProgramGroup missGroup = nullptr;
    OptixProgramGroup hitgroupGroup = nullptr;
};

struct ClusterAccelerationStructureLimits
{
    uint32_t maxTrianglesPerCluster;
    uint32_t maxVerticesPerCluster;
};

struct CUDAMemoryArenaDeleter
{
    void operator()(CUDAMemoryArena *p) const
    {
        if (p)
        {
            p->~CUDAMemoryArena();
            util::AlignedFree(p);
        }
    }
};

struct CUDADevice
{
    CUdevice device;
    CUcontext cudaContext;
    OptixDeviceContext optixDeviceContext;
    size_t totalAllocated;
    size_t bvhTotalAllocated;
    std::unique_ptr<CUDAMemoryArena, CUDAMemoryArenaDeleter> deviceArena;

#if (OPTIX_VERSION >= 90000)
    Array<OptixTraversableHandle> gridClusterTemplateHandles;
#endif

    OptixPrimaryPipelineState optixPrimaryPipeline;

    CUDADevice();
    ~CUDADevice();

    template <typename T>
    DeviceMemoryView<T> Alloc(size_t size);
    template <typename T>
    void Free(DeviceMemoryView<T> &view);
    bool SupportsGrids() const;

    bool CreateOptixPrimaryPipeline(const std::string &ptx);
    void DestroyOptixPrimaryPipeline();

    CUdeviceptr MemAllocBytes(size_t numBytes);
    void MemFreeBytes(CUdeviceptr ptr);
    void MemcpyToDevice(CUdeviceptr dst, const void *src, size_t numBytes);
    void MemcpyToHost(void *dst, CUdeviceptr src, size_t numBytes);

    void BuildBVH(Scene *scene);
    void CreateGridClusterTemplates();
};

OptixDeviceContext InitializeOptix(CUcontext cudaContext);

OptixTraversableHandle
BuildTriangleGASFromMesh(CUDADevice *cudaDevice, HostMemoryArena &hostArena, Mesh &mesh);
OptixTraversableHandle
BuildClusterGASFromMesh(CUDADevice *cudaDevice, HostMemoryArena &hostArena, const Mesh &mesh);
OptixTraversableHandle
BuildCurveGASFromCurves(CUDADevice *cudaDevice,
                        HostMemoryArena &hostArena,
                        Curves &curves);

YBI_NAMESPACE_END

#endif
