#pragma once

#if defined(WITH_CUDA) && defined(WITH_OPTIX)

#include "device/cuda_assert.h"
#include "device/cuda_device_memory_arena.h"
#include "device/device.h"
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

namespace ybi
{

struct Scene;
struct ScenePool;
struct Mesh;
struct Curves;

struct OptixPrimaryPipelineState
{
    OptixPipeline pipeline = nullptr;
    OptixShaderBindingTable sbt = {};
    CUdeviceptr raygenRecordBuffer = 0;
    CUdeviceptr feedbackRaygenRecordBuffer = 0;
    CUdeviceptr missRecordBuffer = 0;
    CUdeviceptr hitgroupRecordBuffer = 0;
    OptixModule module = nullptr;
    OptixModule curveModule = nullptr;
    OptixProgramGroup raygenGroup = nullptr;
    OptixProgramGroup feedbackRaygenGroup = nullptr;
    OptixProgramGroup missGroup = nullptr;
    OptixProgramGroup hitgroupGroup = nullptr;
};

struct ClusterAccelerationStructureLimits
{
    uint32_t maxTrianglesPerCluster = 0;
    uint32_t maxVerticesPerCluster = 0;
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

struct CUDADevice : Device
{
    CUdevice device;
    CUcontext cudaContext;
    OptixDeviceContext optixDeviceContext;
    size_t totalAllocated;
    size_t bvhTotalAllocated;
    std::unique_ptr<CUDAMemoryArena, CUDAMemoryArenaDeleter> deviceArena;

#if (OPTIX_VERSION >= 90000)
    Array<OptixTraversableHandle> gridClusterTemplateHandles;
    bool supportsClusterAccel = false;
    ClusterAccelerationStructureLimits clusterAccelLimits;
#endif

    OptixPrimaryPipelineState optixPrimaryPipeline;

    CUDADevice();
    ~CUDADevice();

    DeviceKind GetKind() const override;
    DeviceMemoryView<uint8_t> AllocBytes(size_t numBytes) override;
    void FreeBytes(DeviceMemoryView<uint8_t> &view) override;
    void CopyBytesToDevice(DeviceMemoryView<uint8_t> dst,
                           const void *src,
                           size_t numBytes) override;
    void CopyBytesToHost(void *dst,
                         DeviceMemoryView<const uint8_t> src,
                         size_t numBytes) override;
    bool InitializeKernels(const std::string &kernelBlob, std::string *outError) override;
    void DestroyKernels() override;
    void ClearTransientMemory() override;
    bool CreateTexture(const DeviceTextureCreateInfo &info,
                       DeviceTexture *outTexture,
                       std::string *outError) override;
    bool UpdateTextureRegion(const DeviceTexture &texture,
                             uint32_t x,
                             uint32_t y,
                             uint32_t width,
                             uint32_t height,
                             const void *pixels,
                             size_t pixelBytes,
                             std::string *outError) override;
    void DestroyTexture(DeviceTexture &texture) override;
    size_t GetBVHAllocatedBytes() const override;

    template <typename T>
    DeviceMemoryView<T> Alloc(size_t size);
    template <typename T>
    void Free(DeviceMemoryView<T> &view);
    bool SupportsGrids() const;
    bool SupportsClusterAccel() const;
    ClusterAccelerationStructureLimits GetClusterAccelerationStructureLimits() const;
    bool DispatchKernel(RenderKernelId kernel, const DispatchParams &params) override;

    bool CreateOptixPrimaryPipeline(const std::string &ptx);
    void DestroyOptixPrimaryPipeline();

    void BuildBVH(Scene *scene) override;
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

} // namespace ybi

#endif
