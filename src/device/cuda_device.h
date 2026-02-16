#pragma once

#include "device/cuda_assert.h"
#include "device/cuda_device_memory_arena.h"
#include "scene/clusterizer.h"
#include "scene/micropolygon_mesh.h"
#include "scene/scene.h"
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

#if defined(WITH_CUDA) && defined(WITH_OPTIX)

struct Scene;

struct ClusterAccelerationStructureLimits
{
    uint32_t maxTrianglesPerCluster;
    uint32_t maxVerticesPerCluster;
};

struct CUDADevice
{
    CUcontext cudaContext;
    OptixDeviceContext optixDeviceContext;
    size_t totalAllocated;
    size_t bvhTotalAllocated;
    std::unique_ptr<CUDAMemoryArena> deviceArena;

#if (OPTIX_VERSION >= 90000)
    Array<OptixTraversableHandle> gridClusterTemplateHandles;
#endif

    CUDADevice();
    ~CUDADevice() = default;

    template <typename T>
    DeviceMemoryView<T> Alloc(size_t size);
    template <typename T>
    void Free(DeviceMemoryView<T> &view);
    bool SupportsGrids() const;

    void BuildBVH(Scene *scene);
    void CreateGridClusterTemplates();
};

#endif

YBI_NAMESPACE_END
