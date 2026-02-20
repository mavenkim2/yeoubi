#include "device/cuda_device.h"
#include "cuda.h"
#include <cstdlib>

YBI_NAMESPACE_BEGIN

#if defined(WITH_CUDA) && defined(WITH_OPTIX)

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

    void *mem = util::AlignedAlloc(sizeof(CUDAMemoryArena), alignof(CUDAMemoryArena));
    YBI_ASSERT(mem != nullptr);
    deviceArena.reset(new (mem) CUDAMemoryArena());
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
    return true;
#else
    return false;
#endif
}

void CUDADevice::CreateGridClusterTemplates()
{
#if (OPTIX_VERSION >= 90000)
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

#endif

YBI_NAMESPACE_END
