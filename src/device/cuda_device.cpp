#include "device/cuda_device.h"
#include "render/dispatch_types.h"
#include "cuda.h"
#include <cstdio>
#include <cstdlib>

YBI_NAMESPACE_BEGIN

#if defined(WITH_CUDA) && defined(WITH_OPTIX)

namespace
{
static size_t KernelIndex(RenderKernelId kernel)
{
    return static_cast<size_t>(kernel);
}

static bool CUDAPrimaryDiffuseStub(CUDADevice *device, const DispatchParams &params)
{
    return device->LaunchPrimaryKernel(params);
}

static bool CUDAAOStub(CUDADevice *device, const DispatchParams &params)
{
    return device->LaunchPrimaryKernel(params);
}

#if (OPTIX_VERSION >= 90000)
static bool QueryOptixUIntProperty(OptixDeviceContext context,
                                   OptixDeviceProperty property,
                                   const char *propertyName,
                                   unsigned int &valueOut)
{
    const OptixResult result =
        optixDeviceContextGetProperty(context, property, &valueOut, sizeof(valueOut));
    if (result == OPTIX_SUCCESS)
    {
        return true;
    }

    fprintf(stderr,
            "OptiX property query failed: %s -> %s (%s)\n",
            propertyName,
            optixGetErrorName(result),
            optixGetErrorString(result));
    return false;
}
#endif
} // namespace

CUDADevice::CUDADevice() : totalAllocated(0), bvhTotalAllocated(0)
{
    CUresult cuResult = cuInit(0);
    if (cuResult != CUDA_SUCCESS)
    {
        const char *str = nullptr;
        const char *name = nullptr;
        cuGetErrorString(cuResult, &str);
        cuGetErrorName(cuResult, &name);
        fprintf(
            stderr, "CUDA init failed: %s (%s)\n", name ? name : "unknown", str ? str : "unknown");
        std::abort();
    }

    cuResult = cuDeviceGet(&device, 0);
    if (cuResult != CUDA_SUCCESS)
    {
        const char *str = nullptr;
        const char *name = nullptr;
        cuGetErrorString(cuResult, &str);
        cuGetErrorName(cuResult, &name);
        fprintf(stderr,
                "CUDA device get failed: %s (%s)\n",
                name ? name : "unknown",
                str ? str : "unknown");
        std::abort();
    }

    cuResult = cuDevicePrimaryCtxRetain(&cudaContext, device);
    if (cuResult != CUDA_SUCCESS)
    {
        const char *str = nullptr;
        const char *name = nullptr;
        cuGetErrorString(cuResult, &str);
        cuGetErrorName(cuResult, &name);
        fprintf(stderr,
                "CUDA primary context retain failed: %s (%s)\n",
                name ? name : "unknown",
                str ? str : "unknown");
        std::abort();
    }

    optixDeviceContext = InitializeOptix(cudaContext);
    if (!optixDeviceContext)
    {
        fprintf(stderr, "Failed to initialize OptiX device context.\n");
        std::abort();
    }

#if (OPTIX_VERSION >= 90000)
    unsigned int clusterAccelFlags = OPTIX_DEVICE_PROPERTY_CLUSTER_ACCEL_FLAG_NONE;
    if (QueryOptixUIntProperty(optixDeviceContext,
                               OPTIX_DEVICE_PROPERTY_CLUSTER_ACCEL,
                               "OPTIX_DEVICE_PROPERTY_CLUSTER_ACCEL",
                               clusterAccelFlags))
    {
        supportsClusterAccel =
            (clusterAccelFlags & OPTIX_DEVICE_PROPERTY_CLUSTER_ACCEL_FLAG_STANDARD) != 0;
    }

    if (supportsClusterAccel)
    {
        unsigned int maxTriangles = 0;
        unsigned int maxVertices = 0;
        const bool gotTriangles = QueryOptixUIntProperty(optixDeviceContext,
                                                         OPTIX_DEVICE_PROPERTY_LIMIT_MAX_CLUSTER_TRIANGLES,
                                                         "OPTIX_DEVICE_PROPERTY_LIMIT_MAX_CLUSTER_TRIANGLES",
                                                         maxTriangles);
        const bool gotVertices = QueryOptixUIntProperty(optixDeviceContext,
                                                        OPTIX_DEVICE_PROPERTY_LIMIT_MAX_CLUSTER_VERTICES,
                                                        "OPTIX_DEVICE_PROPERTY_LIMIT_MAX_CLUSTER_VERTICES",
                                                        maxVertices);
        if (gotTriangles && gotVertices)
        {
            clusterAccelLimits.maxTrianglesPerCluster = maxTriangles;
            clusterAccelLimits.maxVerticesPerCluster = maxVertices;
        }
        else
        {
            supportsClusterAccel = false;
            clusterAccelLimits = {};
            fprintf(stderr,
                    "OptiX cluster acceleration disabled: missing required cluster limits.\n");
        }
    }
#endif

    void *mem = util::AlignedAlloc(sizeof(CUDAMemoryArena), alignof(CUDAMemoryArena));
    YBI_ASSERT(mem != nullptr);
    deviceArena.reset(new (mem) CUDAMemoryArena());

    RegisterKernel(RenderKernelId::PrimaryDiffuse, CUDAPrimaryDiffuseStub);
    RegisterKernel(RenderKernelId::AO, CUDAAOStub);
}

CUDADevice::~CUDADevice()
{
    deviceArena.reset();
    DestroyOptixPrimaryPipeline();
    CUDA_ASSERT(cuDevicePrimaryCtxRelease(device));

    if (optixDeviceContext)
    {
        OPTIX_ASSERT(optixDeviceContextDestroy(optixDeviceContext));
    }
}

DeviceKind CUDADevice::GetKind() const
{
    return DeviceKind::GPU;
}

template <typename T>
DeviceMemoryView<T> CUDADevice::Alloc(size_t count)
{
    YBI_ASSERT(count != 0);
    size_t size = sizeof(T) * count;
    totalAllocated += size;
    CUdeviceptr ptr;
    CUDA_ASSERT(cuMemAlloc(&ptr, size));

    return {(T *)ptr, count};
}

template <typename T>
void CUDADevice::Free(DeviceMemoryView<T> &view)
{
    if (view.data() == nullptr || view.size() == 0)
        return;
    totalAllocated -= view.numBytes();
    CUDA_ASSERT(cuMemFree(CUdeviceptr(view.data())));
    view = {};
}

bool CUDADevice::SupportsGrids() const
{
#if (OPTIX_VERSION >= 90000)
    return supportsClusterAccel;
#else
    return false;
#endif
}

bool CUDADevice::SupportsClusterAccel() const
{
#if (OPTIX_VERSION >= 90000)
    return supportsClusterAccel;
#else
    return false;
#endif
}

ClusterAccelerationStructureLimits CUDADevice::GetClusterAccelerationStructureLimits() const
{
#if (OPTIX_VERSION >= 90000)
    return clusterAccelLimits;
#else
    return {};
#endif
}

DeviceMemoryView<uint8_t> CUDADevice::AllocBytes(size_t numBytes)
{
    YBI_ASSERT(numBytes != 0);
    totalAllocated += numBytes;
    CUdeviceptr ptr = 0;
    CUDA_ASSERT(cuMemAlloc(&ptr, numBytes));
    return {(uint8_t *)ptr, numBytes};
}

void CUDADevice::FreeBytes(DeviceMemoryView<uint8_t> &view)
{
    if (view.data() == nullptr || view.size() == 0)
    {
        return;
    }
    totalAllocated -= view.numBytes();
    CUDA_ASSERT(cuMemFree(CUdeviceptr(view.data())));
    view = {};
}

void CUDADevice::CopyBytesToDevice(DeviceMemoryView<uint8_t> dst,
                                   const void *src,
                                   size_t numBytes)
{
    YBI_ASSERT(dst.data());
    YBI_ASSERT(src);
    YBI_ASSERT(numBytes <= dst.numBytes());
    CUDA_ASSERT(cuMemcpyHtoD(CUdeviceptr(dst.data()), src, numBytes));
}

void CUDADevice::CopyBytesToHost(void *dst,
                                 DeviceMemoryView<const uint8_t> src,
                                 size_t numBytes)
{
    YBI_ASSERT(dst);
    YBI_ASSERT(src.data());
    YBI_ASSERT(numBytes <= src.numBytes());
    CUDA_ASSERT(cuMemcpyDtoH(dst, CUdeviceptr(src.data()), numBytes));
}

size_t CUDADevice::GetBVHAllocatedBytes() const
{
    return bvhTotalAllocated;
}

void CUDADevice::CreateGridClusterTemplates()
{
#if (OPTIX_VERSION >= 90000)
    if (!SupportsClusterAccel())
    {
        fprintf(stderr, "OptiX cluster acceleration not supported by this device/context.\n");
        return;
    }

    const int minDim = 2;
    const int maxDim = 8;
    const int widthDim = maxDim - minDim + 1;
    const int numGrids = widthDim * widthDim;

    unsigned int maxEdges;
    OPTIX_ASSERT(
        optixDeviceContextGetProperty(optixDeviceContext,
                                      OPTIX_DEVICE_PROPERTY_LIMIT_MAX_STRUCTURED_GRID_RESOLUTION,
                                      &maxEdges,
                                      sizeof(unsigned int)));
    YBI_ASSERT(maxDim <= (int)maxEdges);

    OptixClusterAccelBuildInput input = {};
    input.type = OPTIX_CLUSTER_ACCEL_BUILD_TYPE_TEMPLATES_FROM_GRIDS;
    input.grids.flags = OPTIX_CLUSTER_ACCEL_BUILD_FLAG_PREFER_FAST_TRACE;
    input.grids.maxArgCount = (unsigned int)numGrids;
    input.grids.vertexFormat = OPTIX_VERTEX_FORMAT_FLOAT3;
    input.grids.maxSbtIndexValue = 0;
    input.grids.maxWidth = (unsigned int)maxDim;
    input.grids.maxHeight = (unsigned int)maxDim;

    OptixClusterAccelBuildInputGridsArgs hostGridArgs[49];
    int idx = 0;
    for (int w = minDim; w <= maxDim; w++)
    {
        for (int h = minDim; h <= maxDim; h++)
        {
            OptixClusterAccelBuildInputGridsArgs &arg = hostGridArgs[idx];
            arg = {};
            arg.baseClusterId = 0;
            arg.clusterFlags = 0;
            arg.positionTruncateBitCount = 0;
            arg.dimensions[0] = w;
            arg.dimensions[1] = h;
            idx++;
        }
    }

    DeviceMemoryView<OptixClusterAccelBuildInputGridsArgs> args =
        deviceArena->PushArray<OptixClusterAccelBuildInputGridsArgs>((size_t)numGrids);
    DeviceMemoryView<uint32_t> argsCount = deviceArena->PushArray<uint32_t>(1);
    const uint32_t argsCountVal = (uint32_t)numGrids;
    CUDA_ASSERT(cuMemcpyHtoD(CUdeviceptr(args.data()),
                             hostGridArgs,
                             (size_t)numGrids * sizeof(OptixClusterAccelBuildInputGridsArgs)));
    CUDA_ASSERT(cuMemcpyHtoD(CUdeviceptr(argsCount.data()), &argsCountVal, sizeof(uint32_t)));

    OptixAccelBufferSizes bufferSizes = {};
    OPTIX_ASSERT(
        optixClusterAccelComputeMemoryUsage(optixDeviceContext,
                                            OPTIX_CLUSTER_ACCEL_BUILD_MODE_IMPLICIT_DESTINATIONS,
                                            &input,
                                            &bufferSizes));

    YBI_ASSERT(bufferSizes.tempSizeInBytes);
    DeviceMemoryView<uint8_t> temp = deviceArena->PushArray<uint8_t>(
        bufferSizes.tempSizeInBytes, OPTIX_ACCEL_BUFFER_BYTE_ALIGNMENT);

    YBI_ASSERT(bufferSizes.outputSizeInBytes);
    DeviceMemoryView<uint8_t> output = Alloc<uint8_t>(bufferSizes.outputSizeInBytes);

    DeviceMemoryView<OptixTraversableHandle> deviceOutputHandles =
        deviceArena->PushArray<OptixTraversableHandle>((size_t)numGrids);

    OptixClusterAccelBuildModeDesc desc = {};
    desc.mode = OPTIX_CLUSTER_ACCEL_BUILD_MODE_IMPLICIT_DESTINATIONS;
    desc.implicitDest.tempBuffer = CUdeviceptr(temp.data());
    desc.implicitDest.tempBufferSizeInBytes = bufferSizes.tempSizeInBytes;
    desc.implicitDest.outputBuffer = CUdeviceptr(output.data());
    desc.implicitDest.outputBufferSizeInBytes = bufferSizes.outputSizeInBytes;
    desc.implicitDest.outputHandlesBuffer = CUdeviceptr(deviceOutputHandles.data());
    desc.implicitDest.outputHandlesStrideInBytes = sizeof(CUdeviceptr);

    OPTIX_ASSERT(optixClusterAccelBuild(optixDeviceContext,
                                        0,
                                        &desc,
                                        &input,
                                        CUdeviceptr(args.data()),
                                        CUdeviceptr(argsCount.data()),
                                        sizeof(OptixClusterAccelBuildInputGridsArgs)));

    CUDA_ASSERT(cuStreamSynchronize(0));

    gridClusterTemplateHandles.Resize((size_t)numGrids);
    CUDA_ASSERT(cuMemcpyDtoH(gridClusterTemplateHandles.data(),
                             CUdeviceptr(deviceOutputHandles.data()),
                             (size_t)numGrids * sizeof(OptixTraversableHandle)));

    deviceArena->Clear();
#endif
}

bool CUDADevice::DispatchKernel(RenderKernelId kernel, const DispatchParams &params)
{
    const size_t idx = KernelIndex(kernel);
    if (idx >= kernels.size() || kernels[idx] == nullptr)
    {
        fprintf(stderr, "CUDADevice::DispatchKernel missing kernel %zu.\n", idx);
        return false;
    }
    return kernels[idx](this, params);
}

void CUDADevice::RegisterKernel(RenderKernelId kernel, KernelFn fn)
{
    const size_t idx = KernelIndex(kernel);
    YBI_ASSERT(idx < kernels.size());
    kernels[idx] = fn;
}

bool CUDADevice::LaunchPrimaryKernel(const DispatchParams &params)
{
    if (!optixPrimaryPipeline.pipeline)
    {
        fprintf(stderr, "LaunchPrimaryKernel failed: pipeline not initialized.\n");
        return false;
    }
    if (params.launchParamsDevice == 0 || params.launchParamsSize == 0)
    {
        fprintf(stderr, "LaunchPrimaryKernel failed: launch params missing.\n");
        return false;
    }
    if (params.width == 0 || params.height == 0)
    {
        fprintf(stderr, "LaunchPrimaryKernel failed: invalid launch dimensions.\n");
        return false;
    }

    OptixShaderBindingTable sbt = optixPrimaryPipeline.sbt;
    sbt.raygenRecord = optixPrimaryPipeline.raygenRecordBuffer;
    OPTIX_ASSERT(optixLaunch(optixPrimaryPipeline.pipeline,
                             0,
                             static_cast<CUdeviceptr>(params.launchParamsDevice),
                             static_cast<size_t>(params.launchParamsSize),
                             &sbt,
                             params.width,
                             params.height,
                             1));
    CUDA_ASSERT(cuStreamSynchronize(0));
    return true;
}

#endif

YBI_NAMESPACE_END
