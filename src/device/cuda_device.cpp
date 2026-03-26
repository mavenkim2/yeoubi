#include "device/cuda_device.h"
#include "render/dispatch_types.h"
#include <cstdio>
#include <cstdlib>

namespace ybi
{

#if defined(WITH_CUDA) && defined(WITH_OPTIX)

namespace
{
static bool CheckCudaRuntime(cudaError_t result, std::string *outError, const char *callName)
{
    if (result == cudaSuccess)
    {
        return true;
    }
    if (outError)
    {
        *outError = std::string(callName) + " failed: " + cudaGetErrorString(result);
    }
    return false;
}

static cudaTextureAddressMode ToCudaAddressMode(DeviceTextureWrapMode mode)
{
    switch (mode)
    {
        case DeviceTextureWrapMode::Clamp:
            return cudaAddressModeClamp;
        case DeviceTextureWrapMode::Mirror:
            return cudaAddressModeMirror;
        case DeviceTextureWrapMode::Black:
            return cudaAddressModeBorder;
        case DeviceTextureWrapMode::Repeat:
        case DeviceTextureWrapMode::UseMetadata:
        case DeviceTextureWrapMode::Unknown:
        default:
            return cudaAddressModeWrap;
    }
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
}

CUDADevice::~CUDADevice()
{
    deviceArena.reset();
    DestroyKernels();
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

bool CUDADevice::InitializeKernels(const std::string &kernelBlob, std::string *outError)
{
    if (CreateOptixPrimaryPipeline(kernelBlob))
    {
        return true;
    }
    if (outError)
    {
        *outError = "CreateOptixPrimaryPipeline failed.";
    }
    return false;
}

void CUDADevice::DestroyKernels()
{
    DestroyOptixPrimaryPipeline();
}

void CUDADevice::ClearTransientMemory()
{
    if (deviceArena)
    {
        deviceArena->Clear();
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

bool CUDADevice::CreateTexture(const DeviceTextureCreateInfo &info,
                               DeviceTexture *outTexture,
                               std::string *outError)
{
    if (!outTexture)
    {
        if (outError)
        {
            *outError = "CreateTexture: outTexture is null";
        }
        return false;
    }
    *outTexture = {};

    const size_t pixelBytes = TextureFormatPixelBytes(info.format);
    if (pixelBytes == 0u)
    {
        if (outError)
        {
            *outError = "CreateTexture: unsupported format";
        }
        return false;
    }
    if (info.pixels == nullptr || info.width == 0 || info.height == 0)
    {
        if (outError)
        {
            *outError = "CreateTexture: invalid texture input";
        }
        return false;
    }

    const size_t rowBytes = size_t(info.width) * pixelBytes;
    const size_t expectedBytes = rowBytes * size_t(info.height);
    if (info.pixelBytes < expectedBytes)
    {
        if (outError)
        {
            *outError = "CreateTexture: pixelBytes too small";
        }
        return false;
    }

    cudaChannelFormatDesc channelDesc = {};
    cudaTextureReadMode readMode = cudaReadModeElementType;
    switch (info.format)
    {
        case TextureFormat::RGBA8_UNORM:
            channelDesc = cudaCreateChannelDesc(8, 8, 8, 8, cudaChannelFormatKindUnsigned);
            readMode = cudaReadModeNormalizedFloat;
            break;
        case TextureFormat::R16_FLOAT:
            channelDesc = cudaCreateChannelDesc(16, 0, 0, 0, cudaChannelFormatKindFloat);
            readMode = cudaReadModeElementType;
            break;
        case TextureFormat::RG16_FLOAT:
            channelDesc = cudaCreateChannelDesc(16, 16, 0, 0, cudaChannelFormatKindFloat);
            readMode = cudaReadModeElementType;
            break;
        case TextureFormat::RGBA16_FLOAT:
            channelDesc = cudaCreateChannelDesc(16, 16, 16, 16, cudaChannelFormatKindFloat);
            readMode = cudaReadModeElementType;
            break;
        case TextureFormat::R32_FLOAT:
            channelDesc = cudaCreateChannelDesc(32, 0, 0, 0, cudaChannelFormatKindFloat);
            readMode = cudaReadModeElementType;
            break;
        case TextureFormat::RG32_FLOAT:
            channelDesc = cudaCreateChannelDesc(32, 32, 0, 0, cudaChannelFormatKindFloat);
            readMode = cudaReadModeElementType;
            break;
        case TextureFormat::RGBA32_FLOAT:
            channelDesc = cudaCreateChannelDesc(32, 32, 32, 32, cudaChannelFormatKindFloat);
            readMode = cudaReadModeElementType;
            break;
        default:
            if (outError)
            {
                *outError = "CreateTexture: unsupported format";
            }
            return false;
    }
    cudaArray_t array = nullptr;
    if (!CheckCudaRuntime(cudaMallocArray(&array, &channelDesc, info.width, info.height),
                          outError,
                          "cudaMallocArray"))
    {
        return false;
    }

    if (!CheckCudaRuntime(cudaMemcpy2DToArray(array,
                                              0,
                                              0,
                                              info.pixels,
                                              rowBytes,
                                              rowBytes,
                                              info.height,
                                              cudaMemcpyHostToDevice),
                          outError,
                          "cudaMemcpy2DToArray"))
    {
        cudaFreeArray(array);
        return false;
    }

    cudaResourceDesc resourceDesc = {};
    resourceDesc.resType = cudaResourceTypeArray;
    resourceDesc.res.array.array = array;

    cudaTextureDesc textureDesc = {};
    textureDesc.addressMode[0] = ToCudaAddressMode(info.wrapS);
    textureDesc.addressMode[1] = ToCudaAddressMode(info.wrapT);
    textureDesc.filterMode =
        (info.filter == DeviceTextureFilterMode::Linear) ? cudaFilterModeLinear
                                                         : cudaFilterModePoint;
    textureDesc.readMode = readMode;
    textureDesc.normalizedCoords = 1;

    cudaTextureObject_t textureObject = 0;
    if (!CheckCudaRuntime(
            cudaCreateTextureObject(&textureObject, &resourceDesc, &textureDesc, nullptr),
            outError,
            "cudaCreateTextureObject"))
    {
        cudaFreeArray(array);
        return false;
    }

    outTexture->handle = static_cast<uint64_t>(textureObject);
    outTexture->allocation = reinterpret_cast<uint64_t>(array);
    outTexture->width = info.width;
    outTexture->height = info.height;
    outTexture->wrapS = info.wrapS;
    outTexture->wrapT = info.wrapT;
    outTexture->filter = info.filter;
    outTexture->format = info.format;
    outTexture->valid = true;
    return true;
}

bool CUDADevice::UpdateTextureRegion(const DeviceTexture &texture,
                                     uint32_t x,
                                     uint32_t y,
                                     uint32_t width,
                                     uint32_t height,
                                     const void *pixels,
                                     size_t pixelBytes,
                                     std::string *outError)
{
    if (!texture.valid || texture.allocation == 0 || pixels == nullptr)
    {
        if (outError)
        {
            *outError = "UpdateTextureRegion: invalid texture input";
        }
        return false;
    }

    const size_t texelBytes = TextureFormatPixelBytes(texture.format);
    if (texelBytes == 0u || width == 0u || height == 0u || x + width > texture.width ||
        y + height > texture.height)
    {
        if (outError)
        {
            *outError = "UpdateTextureRegion: invalid region";
        }
        return false;
    }

    const size_t rowBytes = static_cast<size_t>(width) * texelBytes;
    const size_t expectedBytes = rowBytes * static_cast<size_t>(height);
    if (pixelBytes < expectedBytes)
    {
        if (outError)
        {
            *outError = "UpdateTextureRegion: pixelBytes too small";
        }
        return false;
    }

    if (!CheckCudaRuntime(cudaMemcpy2DToArray(reinterpret_cast<cudaArray_t>(texture.allocation),
                                              static_cast<size_t>(x) * texelBytes,
                                              y,
                                              pixels,
                                              rowBytes,
                                              rowBytes,
                                              height,
                                              cudaMemcpyHostToDevice),
                          outError,
                          "cudaMemcpy2DToArray"))
    {
        return false;
    }
    return true;
}

void CUDADevice::DestroyTexture(DeviceTexture &texture)
{
    if (texture.handle != 0)
    {
        cudaDestroyTextureObject(static_cast<cudaTextureObject_t>(texture.handle));
    }
    if (texture.allocation != 0)
    {
        cudaFreeArray(reinterpret_cast<cudaArray_t>(texture.allocation));
    }
    texture = {};
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
    if (!optixPrimaryPipeline.pipeline)
    {
        fprintf(stderr, "DispatchKernel failed: pipeline not initialized.\n");
        return false;
    }
    if (params.launchParamsDevice == 0 || params.launchParamsSize == 0)
    {
        fprintf(stderr, "DispatchKernel failed: launch params missing.\n");
        return false;
    }
    if (params.width == 0 || params.height == 0)
    {
        fprintf(stderr, "DispatchKernel failed: invalid launch dimensions.\n");
        return false;
    }

    OptixPipeline optixPipeline = nullptr;
    OptixShaderBindingTable sbt = optixPrimaryPipeline.sbt;

    switch (kernel)
    {
        case RenderKernelId::PrimaryDiffuse:
        case RenderKernelId::AO:
        case RenderKernelId::PathTrace:
            optixPipeline = optixPrimaryPipeline.pipeline;
            sbt.raygenRecord = optixPrimaryPipeline.raygenRecordBuffer;
            break;
        default:
            fprintf(stderr, "DispatchKernel failed: unsupported kernel id.\n");
            return false;
    }

    YBI_ASSERT(optixPipeline != nullptr);
    OPTIX_ASSERT(optixLaunch(optixPipeline,
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

} // namespace ybi
