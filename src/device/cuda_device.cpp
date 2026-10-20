#include "device/cuda_device.h"
#include "core/arena.h"
#include "core/assert.h"
#include "core/bit.h"
#include "core/math.h"
#include "core/thread_context.h"
#include "device/cuda_arena.h"
#include "device/cuda_assert.h"
#include "device/cuda_device.h"
#include "device/device_texture.h"
#include "display_device/display_device.h"
#include "kernel/backend/optix/visibility_mask.h"
#include "kernel/launch/wavefront_params.h"
#include "scene/clusterizer.h"
#include "scene/scene.h"
#include <cstddef>
#include <cstdio>
#include <cstring>
#include <cuda.h>
#include <cuda_profiler_api.h>
#include <cuda_runtime_api.h>
#include <driver_types.h>
#include <fstream>
#include <ios>
#include <istream>
#include <iterator>
#include <memory>
#include <nvtx3/nvtx3.hpp>
#include <nvrtc.h>
#include <optix_function_table_definition.h>
#include <optix_stack_size.h>
#include <optix_stubs.h>
#include <optix_types.h>
#include <string>
#include <vector>

namespace ybi
{

// Helpers

static void PrintLog(const char *log, size_t logSize)
{
    if (log && logSize)
    {
        fprintf(stderr, "%.*s\n", int(logSize), log);
    }
}

#if defined(WITH_CUDA) && defined(WITH_OPTIX)
static OptixResult
CreateOptixModuleCompat(OptixDeviceContext context,
                        const OptixModuleCompileOptions *moduleCompileOptions,
                        const OptixPipelineCompileOptions *pipelineCompileOptions,
                        const char *ptx,
                        size_t ptxSize,
                        char *log,
                        size_t *logSize,
                        OptixModule *moduleOut)
{
#if (OPTIX_VERSION >= 80000)
    return optixModuleCreate(context,
                             moduleCompileOptions,
                             pipelineCompileOptions,
                             ptx,
                             ptxSize,
                             log,
                             logSize,
                             moduleOut);
#else
    return optixModuleCreateFromPTX(context,
                                    moduleCompileOptions,
                                    pipelineCompileOptions,
                                    ptx,
                                    ptxSize,
                                    log,
                                    logSize,
                                    moduleOut);
#endif
}

static OptixResult AccumulateStackSizesCompat(OptixProgramGroup programGroup,
                                              OptixStackSizes *stackSizes,
                                              OptixPipeline pipeline)
{
#if (OPTIX_VERSION >= 80000)
    return optixUtilAccumulateStackSizes(programGroup, stackSizes, pipeline);
#else
    (void)pipeline;
    return optixUtilAccumulateStackSizes(programGroup, stackSizes);
#endif
}

static bool CheckOptixResult(OptixResult result,
                             const char *operation,
                             const char *log = nullptr,
                             size_t logSize = 0)
{
    if (result == OPTIX_SUCCESS)
    {
        return true;
    }

    fprintf(stderr,
            "OptiX Error (%s): %s while calling %s\n",
            optixGetErrorName(result),
            optixGetErrorString(result),
            operation);
    PrintLog(log, logSize);
    return false;
}

static bool CheckCudaResult(CUresult result, const char *operation)
{
    if (result == CUDA_SUCCESS)
        return true;

    const char *name = nullptr;
    const char *message = nullptr;
    cuGetErrorName(result, &name);
    cuGetErrorString(result, &message);
    fprintf(stderr, "CUDA Error (%s): %s while calling %s\n", name, message, operation);
    return false;
}

#endif

void BuildMeshCLAS(CUDADevice *cudaDevice,
                   DataContainer *container,
                   const Mesh &mesh,
                   const MeshletClustersResult &clusterResult,
                   OptixTraversableHandle &outputHandleOut);

// Platform functions
#if defined(WITH_CUDA)

CUDADeviceKernel::CUDADeviceKernel(CUcontext context, CUmodule module, CUfunction function)
    : context(context), module(module), function(function)
{
    CUDA_ASSERT(cuFuncGetName(&name, function));
}

CUDADeviceKernel::~CUDADeviceKernel()
{
    if (!module)
    {
        return;
    }

    CUDA_ASSERT(cuCtxPushCurrent(context));
    CUDA_ASSERT(cuModuleUnload(module));
    CUDA_ASSERT(cuCtxPopCurrent(nullptr));

    module = nullptr;
    function = nullptr;
}

CUDADeviceStream::CUDADeviceStream(CUcontext context) : context(context)
{
    CUDA_ASSERT(cuCtxPushCurrent(context));
    CUDA_ASSERT(cuStreamCreate(&stream, CU_STREAM_NON_BLOCKING));
    CUDA_ASSERT(cuCtxPopCurrent(nullptr));
}

CUDADeviceStream::~CUDADeviceStream()
{
    if (!stream)
    {
        return;
    }

    CUDA_ASSERT(cuCtxPushCurrent(context));
    CUDA_ASSERT(cuStreamDestroy(stream));
    CUDA_ASSERT(cuCtxPopCurrent(nullptr));

    stream = nullptr;
}

CUDADeviceEvent::CUDADeviceEvent(CUcontext context) : context(context)
{
    CUDA_ASSERT(cuCtxPushCurrent(context));
    CUDA_ASSERT(cuEventCreate(&event, CU_EVENT_DISABLE_TIMING));
    CUDA_ASSERT(cuCtxPopCurrent(nullptr));
}

CUDADeviceEvent::~CUDADeviceEvent()
{
    if (!event)
    {
        return;
    }

    CUDA_ASSERT(cuCtxPushCurrent(context));
    CUDA_ASSERT(cuEventDestroy(event));
    CUDA_ASSERT(cuCtxPopCurrent(nullptr));

    event = nullptr;
}

CUDADevice::ContextScope::ContextScope(CUDADevice *device)
{
    CUDA_ASSERT(cuCtxPushCurrent(device->cudaContext));
}

CUDADevice::ContextScope::~ContextScope()
{
    CUDA_ASSERT(cuCtxPopCurrent(0));
}

CUDADevice::CUDADevice(bool createRendererPipeline)
    : device(0), cudaContext(nullptr), optixDeviceContext(nullptr), totalAllocated(0),
      bvhTotalAllocated(0)
{
    CUDA_ASSERT(cuInit(0));
    CUDA_ASSERT(cuDeviceGet(&device, 0));
    CUDA_ASSERT(cuDevicePrimaryCtxRetain(&cudaContext, device));

    InitializeOptix(createRendererPipeline);

    DeviceKernelDesc cameraDesc;
    cameraDesc.filename = YBI_GENERATE_CAMERA_RAY_KERNEL_PATH;
    cameraDesc.entryName = "GenerateCameraRays";
    generateCameraRays = CreateDeviceKernel(cameraDesc);
}

CUDADevice::~CUDADevice()
{
    {
        ContextScope contextScope(this);
        generateCameraRays.reset();
        DestroyOptixPipeline();
        if (optixDeviceContext)
        {
            OPTIX_ASSERT(optixDeviceContextDestroy(optixDeviceContext));
            optixDeviceContext = nullptr;
        }
    }
    CUDA_ASSERT(cuDevicePrimaryCtxRelease(device));
    cudaContext = nullptr;
}

void CUDADevice::WaitIdle()
{
    ContextScope contextScope(this);
    CUDA_ASSERT(cuCtxSynchronize());
}

DeviceBuffer<std::byte> CUDADevice::AllocBytes(size_t numBytes)
{
    YBI_ASSERT(numBytes != 0);
    size_t size = numBytes;
    // totalAllocated += size;
    CUdeviceptr ptr;
    CUDA_ASSERT(cuMemAlloc(&ptr, size));

    return DeviceBuffer<std::byte>(uintptr_t(ptr), numBytes);
}

void CUDADevice::FreeBytes(DeviceBuffer<std::byte> &buffer)
{
    if (buffer.Data() == nullptr || buffer.SizeInBytes() == 0)
        return;

    // totalAllocated -= view.numBytes();
    CUDA_ASSERT(cuMemFree(CUdeviceptr(buffer.Data())));

    buffer = DeviceBuffer<std::byte>(std::nullptr_t());
}

DeviceBuffer<std::byte> CUDADevice::AllocHostBytes(size_t numBytes)
{
    YBI_ASSERT(numBytes != 0);

    void *ptr = nullptr;
    CUDA_ASSERT(cuMemAllocHost(&ptr, numBytes));
    return DeviceBuffer<std::byte>(static_cast<std::byte *>(ptr), numBytes);
}

void CUDADevice::FreeHostBytes(DeviceBuffer<std::byte> &buffer)
{
    if (buffer.Data() == nullptr || buffer.SizeInBytes() == 0)
        return;

    CUDA_ASSERT(cuMemFreeHost(buffer.Data()));
    buffer = DeviceBuffer<std::byte>(std::nullptr_t());
}

DeviceExternalBuffer CUDADevice::ImportExternalBuffer(const DisplayBuffer &buffer)
{
    YBI_ASSERT(buffer.size > 0);
    YBI_ASSERT(buffer.memory.size > 0);

    CUDA_EXTERNAL_MEMORY_HANDLE_DESC desc = {};
    switch (buffer.memory.type)
    {
        case ExternalMemoryHandleType::OpaqueWin32:
        {
            YBI_ASSERT(buffer.memory.handle);
            desc.type = CU_EXTERNAL_MEMORY_HANDLE_TYPE_OPAQUE_WIN32;
            desc.handle.win32.handle = buffer.memory.handle;
            break;
        }
        default:
        {
            YBI_ASSERT(false);
            break;
        }
    }

    desc.size = buffer.memory.size;

    CUexternalMemory externalMemory;
    CUDA_ASSERT(cuImportExternalMemory(&externalMemory, &desc));

    CUDA_EXTERNAL_MEMORY_BUFFER_DESC bufferDesc = {};
    bufferDesc.offset = 0;
    bufferDesc.size = buffer.size;

    CUdeviceptr cudaPtr = 0;
    CUDA_ASSERT(cuExternalMemoryGetMappedBuffer(&cudaPtr, externalMemory, &bufferDesc));

    // importedExternalBuffers.push_back({cudaPtr, buffer.size, externalMemory});
    DeviceExternalBuffer result;
    result.handle = reinterpret_cast<void *>(externalMemory);
    result.buffer = DeviceBuffer<std::byte>(uintptr_t(cudaPtr), buffer.size);
    result.device = this;
    return result;
}

void CUDADevice::FreeExternalBuffer(DeviceExternalBuffer &buffer)
{
    if (buffer.buffer)
    {
        CUDA_ASSERT(cuMemFree(reinterpret_cast<CUdeviceptr>(buffer.buffer.Data())));
        buffer.buffer = DeviceBuffer<std::byte>(std::nullptr_t());
    }

    if (buffer.handle)
    {
        CUDA_ASSERT(cuDestroyExternalMemory(reinterpret_cast<CUexternalMemory>(buffer.handle)));
        buffer.handle = nullptr;
    }
}

DeviceExternalSemaphore CUDADevice::ImportExternalSemaphore(const DisplaySemaphore &semaphore)
{
    CUDA_EXTERNAL_SEMAPHORE_HANDLE_DESC desc = {};
    switch (semaphore.handle.type)
    {
        case ExternalMemoryHandleType::OpaqueWin32:
        {
            YBI_ASSERT(semaphore.handle.handle);
            desc.type = CU_EXTERNAL_SEMAPHORE_HANDLE_TYPE_OPAQUE_WIN32;
            desc.handle.win32.handle = semaphore.handle.handle;
            break;
        }
        default:
        {
            YBI_ASSERT(false);
            break;
        }
    }

    CUexternalSemaphore externalSemaphore = nullptr;
    CUDA_ASSERT(cuImportExternalSemaphore(&externalSemaphore, &desc));

    DeviceExternalSemaphore result;
    result.handle = reinterpret_cast<void *>(externalSemaphore);
    result.device = this;
    return result;
}

void CUDADevice::FreeExternalSemaphore(DeviceExternalSemaphore &semaphore)
{
    if (!semaphore.handle)
    {
        return;
    }

    CUDA_ASSERT(
        cuDestroyExternalSemaphore(reinterpret_cast<CUexternalSemaphore>(semaphore.handle)));
    semaphore.handle = nullptr;
}

void CUDADevice::SignalExternalSemaphore(DeviceExternalSemaphore &semaphore)
{
    YBI_ASSERT(semaphore.handle);

    CUexternalSemaphore externalSemaphore =
        reinterpret_cast<CUexternalSemaphore>(semaphore.handle);
    CUDA_EXTERNAL_SEMAPHORE_SIGNAL_PARAMS signalParams = {};
    CUDA_ASSERT(cuSignalExternalSemaphoresAsync(&externalSemaphore, &signalParams, 1, 0));
}

std::unique_ptr<DeviceArena> CUDADevice::CreateArena()
{
    return std::make_unique<CUDAArena>();
}

void CUDADevice::CopyToHost(void *dst, void *ptr, size_t size)
{
    CUDA_ASSERT(cuMemcpyDtoH(dst, reinterpret_cast<CUdeviceptr>(ptr), size));
}

void CUDADevice::CopyToHostAsync(void *dst, void *ptr, size_t size)
{
    CUDA_ASSERT(cuMemcpyDtoHAsync(dst, reinterpret_cast<CUdeviceptr>(ptr), size, 0));
}

void CUDADevice::CopyToDevice(void *dst, void *ptr, size_t size)
{
    CUDA_ASSERT(cuMemcpyHtoD(reinterpret_cast<CUdeviceptr>(dst), ptr, size));
}

void CUDADevice::CopyToDeviceAsync(void *dst, void *ptr, size_t size)
{
    CUDA_ASSERT(cuMemcpyHtoDAsync(reinterpret_cast<CUdeviceptr>(dst), ptr, size, 0));
}

void CUDADevice::CopyDeviceToDeviceAsync(void *dst, void *ptr, size_t size)
{
    ContextScope scope(this);
    CUDA_ASSERT(cuMemcpyDtoDAsync(reinterpret_cast<CUdeviceptr>(dst),
                                reinterpret_cast<CUdeviceptr>(ptr), size, 0));
}

void CUDADevice::Memzero(void *dst, size_t size)
{
    CUDA_ASSERT(cuMemsetD8(reinterpret_cast<CUdeviceptr>(dst), 0, size));
}

void CUDADevice::MemzeroAsync(void *dst, size_t size)
{
    CUDA_ASSERT(cuMemsetD8Async(reinterpret_cast<CUdeviceptr>(dst), 0, size, 0));
}

std::unique_ptr<DeviceStream> CUDADevice::CreateDeviceStream()
{
    return std::make_unique<CUDADeviceStream>(cudaContext);
}

std::shared_ptr<DeviceEvent> CUDADevice::CreateDeviceEvent()
{
    return std::make_shared<CUDADeviceEvent>(cudaContext);
}

void CUDADevice::SyncDeviceStream(DeviceStream *stream)
{
    ContextScope contextScope(this);

    CUstream cudaStream = nullptr;
    if (stream)
    {
        CUDADeviceStream &deviceStream = static_cast<CUDADeviceStream &>(*stream);
        YBI_ASSERT(deviceStream.stream);
        cudaStream = deviceStream.stream;
    }
    CUDA_ASSERT(cuStreamSynchronize(cudaStream));
}

void CUDADevice::RecordDeviceEvent(DeviceStream *stream, DeviceEvent &event)
{
    ContextScope contextScope(this);

    CUDADeviceEvent &cudaEvent = static_cast<CUDADeviceEvent &>(event);
    YBI_ASSERT(cudaEvent.event);

    CUstream cudaStream = nullptr;
    if (stream)
    {
        CUDADeviceStream &deviceStream = static_cast<CUDADeviceStream &>(*stream);
        YBI_ASSERT(deviceStream.stream);
        cudaStream = deviceStream.stream;
    }

    CUDA_ASSERT(cuEventRecord(cudaEvent.event, cudaStream));
}

void CUDADevice::WaitDeviceEvent(DeviceStream &stream, DeviceEvent &event)
{
    ContextScope contextScope(this);

    CUDADeviceStream &cudaStream = static_cast<CUDADeviceStream &>(stream);
    CUDADeviceEvent &cudaEvent = static_cast<CUDADeviceEvent &>(event);
    YBI_ASSERT(cudaStream.stream);
    YBI_ASSERT(cudaEvent.event);

    CUDA_ASSERT(cuStreamWaitEvent(cudaStream.stream, cudaEvent.event, 0));
}

void CUDADevice::CopyToDevice3DAsync(DeviceStream &stream,
                                     DeviceTexture &texture,
                                     const DeviceTextureCopyDesc &desc)
{
    ContextScope contextScope(this);

    CUDADeviceStream &cudaStream = static_cast<CUDADeviceStream &>(stream);
    YBI_ASSERT(cudaStream.stream);
    YBI_ASSERT(texture.array);
    YBI_ASSERT(desc.src);
    YBI_ASSERT(desc.srcRowPitch >= desc.widthInBytes);
    YBI_ASSERT(desc.widthInBytes > 0);
    YBI_ASSERT(desc.height > 0);
    YBI_ASSERT(desc.depth > 0);

    const uint32_t elementSize = GetTextureFormatSize(texture.info.format);
    const uint32_t textureDepth =
        Max(1u, texture.info.depth ? texture.info.depth : texture.info.layers);
    YBI_ASSERT(elementSize > 0);
    YBI_ASSERT(desc.widthInBytes % elementSize == 0);
    YBI_ASSERT(desc.dstX + (desc.widthInBytes / elementSize) <= texture.info.width);
    YBI_ASSERT(desc.dstY + desc.height <= Max(1u, texture.info.height));
    YBI_ASSERT(desc.dstZ + desc.depth <= textureDepth);
    YBI_ASSERT(desc.srcSlicePitch == 0 || desc.srcSlicePitch >= desc.srcRowPitch * desc.height);
    YBI_ASSERT(desc.srcSlicePitch == 0 || desc.srcSlicePitch % desc.srcRowPitch == 0);

    CUDA_MEMCPY3D copy = {};
    copy.srcMemoryType = CU_MEMORYTYPE_HOST;
    copy.srcHost = desc.src;
    copy.srcPitch = desc.srcRowPitch;
    copy.srcHeight = desc.srcSlicePitch ? desc.srcSlicePitch / desc.srcRowPitch : desc.height;

    copy.dstMemoryType = CU_MEMORYTYPE_ARRAY;
    copy.dstArray = reinterpret_cast<CUarray>(texture.array);
    copy.dstXInBytes = desc.dstX * elementSize;
    copy.dstY = desc.dstY;
    copy.dstZ = desc.dstZ;

    copy.WidthInBytes = desc.widthInBytes;
    copy.Height = desc.height;
    copy.Depth = desc.depth;

    CUDA_ASSERT(cuMemcpy3DAsync(&copy, cudaStream.stream));
}

void CUDADevice::CopyToDevice3DBatchAsync(DeviceStream &stream,
                                          const DeviceTextureBatchCopyDesc *descs,
                                          size_t numDescs)
{
    if (numDescs == 0)
    {
        return;
    }

    ContextScope contextScope(this);

    CUDADeviceStream &cudaStream = static_cast<CUDADeviceStream &>(stream);
    YBI_ASSERT(cudaStream.stream);
    YBI_ASSERT(descs);

    std::vector<CUDA_MEMCPY3D_BATCH_OP> ops(numDescs);
    for (size_t descIndex = 0; descIndex < numDescs; descIndex++)
    {
        const DeviceTextureBatchCopyDesc &batchDesc = descs[descIndex];
        YBI_ASSERT(batchDesc.texture);

        DeviceTexture &texture = *batchDesc.texture;
        const DeviceTextureCopyDesc &copyDesc = batchDesc.copy;

        YBI_ASSERT(texture.array);
        YBI_ASSERT(copyDesc.src);
        YBI_ASSERT(copyDesc.srcRowPitch >= copyDesc.widthInBytes);
        YBI_ASSERT(copyDesc.widthInBytes > 0);
        YBI_ASSERT(copyDesc.height > 0);
        YBI_ASSERT(copyDesc.depth > 0);

        const uint32_t elementSize = GetTextureFormatSize(texture.info.format);
        const uint32_t textureDepth =
            Max(1u, texture.info.depth ? texture.info.depth : texture.info.layers);
        YBI_ASSERT(elementSize > 0);
        YBI_ASSERT(copyDesc.srcRowPitch % elementSize == 0);
        YBI_ASSERT(copyDesc.widthInBytes % elementSize == 0);
        YBI_ASSERT(copyDesc.dstX + (copyDesc.widthInBytes / elementSize) <= texture.info.width);
        YBI_ASSERT(copyDesc.dstY + copyDesc.height <= Max(1u, texture.info.height));
        YBI_ASSERT(copyDesc.dstZ + copyDesc.depth <= textureDepth);
        YBI_ASSERT(copyDesc.srcSlicePitch == 0 ||
                   copyDesc.srcSlicePitch >= copyDesc.srcRowPitch * copyDesc.height);
        YBI_ASSERT(copyDesc.srcSlicePitch == 0 ||
                   copyDesc.srcSlicePitch % copyDesc.srcRowPitch == 0);

        CUDA_MEMCPY3D_BATCH_OP &op = ops[descIndex];
        op = {};

        op.src.type = CU_MEMCPY_OPERAND_TYPE_POINTER;
        op.src.op.ptr.ptr = reinterpret_cast<CUdeviceptr>(copyDesc.src);
        op.src.op.ptr.rowLength = copyDesc.srcRowPitch / elementSize;
        op.src.op.ptr.layerHeight = copyDesc.srcSlicePitch
                                        ? copyDesc.srcSlicePitch / copyDesc.srcRowPitch
                                        : copyDesc.height;

        op.dst.type = CU_MEMCPY_OPERAND_TYPE_ARRAY;
        op.dst.op.array.array = reinterpret_cast<CUarray>(texture.array);
        op.dst.op.array.offset.x = copyDesc.dstX;
        op.dst.op.array.offset.y = copyDesc.dstY;
        op.dst.op.array.offset.z = copyDesc.dstZ;

        op.extent.width = copyDesc.widthInBytes / elementSize;
        op.extent.height = copyDesc.height;
        op.extent.depth = copyDesc.depth;
        op.srcAccessOrder = CU_MEMCPY_SRC_ACCESS_ORDER_STREAM;
    }

    CUDA_ASSERT(cuMemcpy3DBatchAsync(ops.size(), ops.data(), 0, cudaStream.stream));
}

bool CUDADevice::DispatchKernel(KernelType kernelType,
                                uint32_t numWorkItems,
                                uint32_t blockSize,
                                DeviceBuffer<std::byte> kernelParams)
{
    ContextScope contextScope(this);
    if (!kernelParams.Data() || kernelParams.SizeInBytes() < sizeof(WavefrontParams))
        return false;
    if (numWorkItems == 0)
        return true;

    if (kernelType == KERNEL_GENERATE_CAMERA_RAYS)
    {
        const nvtx3::scoped_range range{"GenerateCameraRays"};
        if (blockSize == 0)
            return false;

        CUDADeviceKernel &cameraKernel = static_cast<CUDADeviceKernel &>(*generateCameraRays);
        void *arguments[] = {kernelParams.KernelArg(), &numWorkItems};
        CUresult result = cuLaunchKernel(cameraKernel.function,
                                         CeilDiv(numWorkItems, blockSize),
                                         1,
                                         1,
                                         blockSize,
                                         1,
                                         1,
                                         0,
                                         0,
                                         arguments,
                                         nullptr);
        if (!CheckCudaResult(result, "cuLaunchKernel(GenerateCameraRays)"))
            return false;

        return true;
    }

    if (kernelType < KERNEL_INTERSECT || kernelType > KERNEL_INDIRECT_LIGHTING ||
        !optixPipeline.pipeline)
        return false;
    const uint32_t raygen = kernelType - KERNEL_INTERSECT;
    static constexpr const char *names[] = {
        "Intersect", "ShadeHit", "DirectLighting", "IndirectLighting"};
    const nvtx3::scoped_range range{names[raygen]};
    OptixShaderBindingTable sbt = {};
    sbt.raygenRecord = optixPipeline.sbtRecordsBuffer + sizeof(PipelineSbtRecord) * raygen;
    sbt.missRecordBase =
        optixPipeline.sbtRecordsBuffer + sizeof(PipelineSbtRecord) * PROGRAM_GROUP_MISS;
    sbt.missRecordStrideInBytes = sizeof(PipelineSbtRecord);
    sbt.missRecordCount = 1;
    sbt.hitgroupRecordBase =
        optixPipeline.sbtRecordsBuffer + sizeof(PipelineSbtRecord) * PROGRAM_GROUP_CLOSEST_HIT;
    sbt.hitgroupRecordStrideInBytes = sizeof(PipelineSbtRecord);
    sbt.hitgroupRecordCount = 1;
    if (optixPipeline.callableRecords)
    {
        sbt.callablesRecordBase = optixPipeline.callableRecords;
        sbt.callablesRecordStrideInBytes = sizeof(PipelineSbtRecord);
        sbt.callablesRecordCount =
            static_cast<unsigned int>(optixPipeline.callableProgramGroups.size());
    }

    OptixResult result = optixLaunch(optixPipeline.pipeline,
                                     0,
                                     reinterpret_cast<CUdeviceptr>(kernelParams.Data()),
                                     sizeof(WavefrontParams),
                                     &sbt,
                                     numWorkItems,
                                     1,
                                     1);
    if (!CheckOptixResult(result, "optixLaunch"))
        return false;

    return true;
}

std::unique_ptr<DeviceKernel> CUDADevice::CreateDeviceKernel(const DeviceKernelDesc &desc)
{
    ContextScope contextScope(this);

    YBI_ASSERT(!desc.filename.empty());

    CUmodule module = nullptr;
    CUDA_ASSERT(cuModuleLoad(&module, desc.filename.c_str()));

    CUfunction function = nullptr;
    if (!desc.entryName.empty())
    {
        CUDA_ASSERT(cuModuleGetFunction(&function, module, desc.entryName.c_str()));
    }
    else
    {
        uint32_t numFunctions = 0;
        CUDA_ASSERT(cuModuleGetFunctionCount(&numFunctions, module));
        YBI_ASSERT(numFunctions == 1);

        CUDA_ASSERT(cuModuleEnumerateFunctions(&function, 1, module));
    }

    return std::make_unique<CUDADeviceKernel>(cudaContext, module, function);
}

bool CUDADevice::CompileKernel(const std::string &source, DeviceModule &deviceModule)
{
    nvrtcProgram program;

    nvrtcResult nvrtcStatus =
        nvrtcCreateProgram(&program, source.c_str(), nullptr, 0, nullptr, nullptr);

    YBI_ASSERT(nvrtcStatus == NVRTC_SUCCESS);

    std::vector<std::string> optionStorage = {
        "--std=" YBI_CUDA_STANDARD,
        "--gpu-architecture=" YBI_CUDA_ARCHITECTURE,
        "--device-as-default-execution-space",
        "--relocatable-device-code=true",
        "--use_fast_math",
        "--generate-line-info",
        "-DOPTIX_VERSION=" YBI_OPTIX_VERSION_STRING,
#if defined(YBI_USE_PATH_GUIDING)
        "-DYBI_USE_PATH_GUIDING",
#endif
        "-I" YBI_OPTIX_INCLUDE_PATH,
        "-I" YBI_OPTIX_SDK_INCLUDE_PATH,
        "-I" YBI_CUDA_INCLUDE_PATH,
#if defined(YBI_CUDA_CCCL_INCLUDE_PATH)
        "-I" YBI_CUDA_CCCL_INCLUDE_PATH,
#endif
    };

    std::vector<const char *> options;
    options.reserve(optionStorage.size());
    for (const std::string &option : optionStorage)
        options.push_back(option.c_str());

    nvrtcResult compileResult =
        nvrtcCompileProgram(program, static_cast<int>(options.size()), options.data());

    size_t logSize = 0;
    nvrtcStatus = nvrtcGetProgramLogSize(program, &logSize);

    std::vector<char> log(Max(4096llu, logSize));
    if (logSize)
    {
        nvrtcStatus = nvrtcGetProgramLog(program, log.data());
    }

    if (compileResult != NVRTC_SUCCESS)
    {
        PrintLog(log.data(), logSize);
        nvrtcDestroyProgram(&program);
        return false;
    }

    size_t ptxSize = 0;
    nvrtcGetPTXSize(program, &ptxSize);

    std::vector<char> ptx(ptxSize);
    nvrtcGetPTX(program, ptx.data());

    nvrtcDestroyProgram(&program);

    YBI_ERROR(optixDeviceContext, "Must initialize OptiX device context before compiling kernel");

#if defined(WITH_OPTIX)
    ContextScope contextScope(this);

    OptixModuleCompileOptions moduleCompileOptions = {};
    OptixPipelineCompileOptions pipelineCompileOptions = {};

    GetOptixCompileOptions(moduleCompileOptions, pipelineCompileOptions);

    OptixModule module = nullptr;
    logSize = log.size();
    OptixResult result = CreateOptixModuleCompat(optixDeviceContext,
                                                 &moduleCompileOptions,
                                                 &pipelineCompileOptions,
                                                 ptx.data(),
                                                 ptx.size(),
                                                 log.data(),
                                                 &logSize,
                                                 &module);
    if (!CheckOptixResult(result, "optixModuleCreate", log.data(), logSize))
    {
        return false;
    }

    deviceModule = reinterpret_cast<DeviceModule>(module);
    return true;
#else
#error cuda w/out optix not supported
    return false;
#endif
}

bool CUDADevice::RegisterDirectCallables(DeviceModule deviceModule,
                                         const std::vector<std::string> &entryNames,
                                         std::vector<uint32_t> &callableIndices)
{
#if defined(WITH_OPTIX)
    ContextScope contextScope(this);

    callableIndices.clear();
    if (deviceModule == 0 || entryNames.empty() ||
        !optixPipeline.programGroups[PROGRAM_GROUP_RAYGEN_INTERSECT])
        return false;

    const OptixModule module = reinterpret_cast<OptixModule>(deviceModule);
    const size_t firstCallableIndex = optixPipeline.callableProgramGroups.size();
    OptixProgramGroupOptions programGroupOptions = {};
    char log[4096] = {};

    for (const std::string &entryName : entryNames)
    {
        OptixProgramGroupDesc desc = {};
        desc.kind = OPTIX_PROGRAM_GROUP_KIND_CALLABLES;
        desc.callables.moduleDC = module;
        desc.callables.entryFunctionNameDC = entryName.c_str();

        OptixProgramGroup programGroup = nullptr;
        size_t logSize = sizeof(log);
        OptixResult result = optixProgramGroupCreate(
            optixDeviceContext, &desc, 1, &programGroupOptions, log, &logSize, &programGroup);
        if (!CheckOptixResult(result, "optixProgramGroupCreate(callable)", log, logSize))
        {
            for (size_t index = firstCallableIndex;
                 index < optixPipeline.callableProgramGroups.size();
                 ++index)
            {
                result = optixProgramGroupDestroy(optixPipeline.callableProgramGroups[index]);
                CheckOptixResult(result, "optixProgramGroupDestroy(callable)");
            }
            optixPipeline.callableProgramGroups.resize(firstCallableIndex);
            result = optixModuleDestroy(module);
            CheckOptixResult(result, "optixModuleDestroy(callable)");
            return false;
        }

        optixPipeline.callableProgramGroups.push_back(programGroup);
    }

    optixPipeline.callableModules.push_back(module);

    OptixModuleCompileOptions moduleCompileOptions = {};
    OptixPipelineCompileOptions pipelineCompileOptions = {};
    GetOptixCompileOptions(moduleCompileOptions, pipelineCompileOptions);
    if (!LinkOptixPipeline(pipelineCompileOptions))
        return false;

    callableIndices.reserve(entryNames.size());
    for (uint32_t index = 0; index < static_cast<uint32_t>(entryNames.size()); ++index)
    {
        callableIndices.push_back(static_cast<uint32_t>(firstCallableIndex) + index);
    }
    return true;
#else
    (void)deviceModule;
    (void)entryNames;
    (void)callableIndices;
    return false;
#endif
}

std::string CUDADevice::GetCallableFunctionSignature() const
{
    return "extern \"C\" __device__";
}

std::string CUDADevice::GetDirectCallablePrefix() const
{
    return "__direct_callable__";
}

void CUDADevice::ExecuteKernel(DeviceKernel &kernel,
                               const DeviceKernelLaunchDesc &desc)
{
    ContextScope contextScope(this);

    CUDADeviceKernel &cudaKernel = static_cast<CUDADeviceKernel &>(kernel);
    const nvtx3::scoped_range range{cudaKernel.name};
    YBI_ASSERT(cudaKernel.function);
    YBI_ASSERT(desc.arguments);
    YBI_ASSERT(desc.gridX > 0);
    YBI_ASSERT(desc.gridY > 0);
    YBI_ASSERT(desc.gridZ > 0);
    YBI_ASSERT(desc.blockX > 0);
    YBI_ASSERT(desc.blockY > 0);
    YBI_ASSERT(desc.blockZ > 0);

    CUDA_ASSERT(cuLaunchKernel(cudaKernel.function,
                               desc.gridX,
                               desc.gridY,
                               desc.gridZ,
                               desc.blockX,
                               desc.blockY,
                               desc.blockZ,
                               desc.sharedMemoryBytes,
                               0,
                               desc.arguments,
                               nullptr));
}

struct CUDATextureFormatInfo
{
    CUarray_format format = CU_AD_FORMAT_UNSIGNED_INT8;
    uint32_t numChannels = 0;
};

static CUDATextureFormatInfo GetCUDATextureFormatInfo(TextureFormat format)
{
    CUDATextureFormatInfo result = {};

    switch (format)
    {
        case TextureFormat::R8_UNORM:
        case TextureFormat::RG8_UNORM:
        case TextureFormat::RGBA8_UNORM:
        {
            result.format = CU_AD_FORMAT_UNSIGNED_INT8;
            result.numChannels = GetTextureFormatNumChannels(format);
            break;
        }
        case TextureFormat::R16_UNORM:
        case TextureFormat::RG16_UNORM:
        case TextureFormat::RGBA16_UNORM:
        {
            result.format = CU_AD_FORMAT_UNSIGNED_INT16;
            result.numChannels = GetTextureFormatNumChannels(format);
            break;
        }
        case TextureFormat::R16_FLOAT:
        {
            result.format = CU_AD_FORMAT_HALF;
            result.numChannels = 1;
            break;
        }
        case TextureFormat::RG16_FLOAT:
        {
            result.format = CU_AD_FORMAT_HALF;
            result.numChannels = 2;
            break;
        }
        case TextureFormat::RGBA16_FLOAT:
        {
            result.format = CU_AD_FORMAT_HALF;
            result.numChannels = 4;
            break;
        }
        case TextureFormat::R32_FLOAT:
        {
            result.format = CU_AD_FORMAT_FLOAT;
            result.numChannels = 1;
            break;
        }
        case TextureFormat::RG32_FLOAT:
        {
            result.format = CU_AD_FORMAT_FLOAT;
            result.numChannels = 2;
            break;
        }
        case TextureFormat::RGBA32_FLOAT:
        {
            result.format = CU_AD_FORMAT_FLOAT;
            result.numChannels = 4;
            break;
        }
        default:
        {
            YBI_ASSERT(false);
            break;
        }
    }

    return result;
}

DeviceTexture CUDADevice::CreateTexture(const DeviceTextureInfo &info)
{
    ContextScope contextScope(this);

    YBI_ASSERT(info.width > 0);
    YBI_ASSERT(info.height > 0 || (info.layers == 0 && info.depth == 0));
    YBI_ASSERT(info.layers == 0 || info.depth == 0);

    const CUDATextureFormatInfo formatInfo = GetCUDATextureFormatInfo(info.format);

    CUfilter_mode filterMode;
    switch (info.filterType)
    {
        case FilterType::Point:
        {
            filterMode = CU_TR_FILTER_MODE_POINT;
            break;
        }
        case FilterType::Linear:
        {
            filterMode = CU_TR_FILTER_MODE_LINEAR;
            break;
        }
    }

    CUarray array = nullptr;
    CUDA_ARRAY3D_DESCRIPTOR descriptor = {};
    descriptor.Width = info.width;
    descriptor.Height = info.height;
    descriptor.Depth = info.depth ? info.depth : info.layers;
    descriptor.Format = formatInfo.format;
    descriptor.NumChannels = formatInfo.numChannels;
    descriptor.Flags = info.layers ? CUDA_ARRAY3D_LAYERED : 0;
    // if (info.enableSurfaceWrites)
    // {
    //     descriptor.Flags |= CUDA_ARRAY3D_SURFACE_LDST;
    // }

    CUDA_ASSERT(cuArray3DCreate(&array, &descriptor));

    CUDA_RESOURCE_DESC resourceDesc = {};
    resourceDesc.resType = CU_RESOURCE_TYPE_ARRAY;
    resourceDesc.res.array.hArray = array;

    CUDA_TEXTURE_DESC textureDesc = {};
    for (uint32_t axis = 0; axis < 3; axis++)
    {
        textureDesc.addressMode[axis] = info.addressMode[axis] == TextureAddressMode::Wrap
                                            ? CU_TR_ADDRESS_MODE_WRAP
                                            : CU_TR_ADDRESS_MODE_CLAMP;
    }
    textureDesc.filterMode = filterMode;
    textureDesc.flags = CU_TRSF_NORMALIZED_COORDINATES;
    textureDesc.maxAnisotropy = 1;
    textureDesc.mipmapFilterMode = CU_TR_FILTER_MODE_POINT;
    textureDesc.minMipmapLevelClamp = 0.0f;
    textureDesc.maxMipmapLevelClamp = info.mips ? float(info.mips - 1) : 0;

    CUtexObject textureObject = 0;
    CUDA_ASSERT(cuTexObjectCreate(&textureObject, &resourceDesc, &textureDesc, nullptr));

    // CUsurfObject surfaceObject = 0;
    // if (info.enableSurfaceWrites)
    // {
    //     CUDA_ASSERT(cuSurfObjectCreate(&surfaceObject, &resourceDesc));
    // }

    DeviceTexture result = {};
    result.array = reinterpret_cast<uintptr_t>(array);
    result.textureObject = textureObject;
    result.info = info;
    // result.surfaceObject = surfaceObject;
    return result;
}

void CUDADevice::DestroyTexture(DeviceTexture &texture)
{
    if (texture.textureObject == 0 && texture.array == 0)
        return;

    ContextScope contextScope(this);
    if (texture.textureObject != 0)
    {
        CUDA_ASSERT(cuTexObjectDestroy(static_cast<CUtexObject>(texture.textureObject)));
    }
    if (texture.array != 0)
    {
        CUDA_ASSERT(cuArrayDestroy(reinterpret_cast<CUarray>(texture.array)));
    }

    texture = {};
}

#endif

// BVH
#if defined(WITH_CUDA) && defined(WITH_OPTIX)

void CUDADevice::InitializeOptix(bool createRendererPipeline)
{
    OPTIX_ASSERT(optixInit());

    OptixDeviceContextOptions contextOptions = {};
    contextOptions.logCallbackFunction =
        [](unsigned int level, const char *tag, const char *message, void *cbdata) {
            fprintf(stderr, "[OptiX][%u][%s] %s\n", level, tag, message);
            (void)cbdata;
        };
    contextOptions.logCallbackLevel = 4;
    contextOptions.validationMode = OPTIX_DEVICE_CONTEXT_VALIDATION_MODE_ALL;

    OPTIX_ASSERT(optixDeviceContextCreate(cudaContext, &contextOptions, &optixDeviceContext));

    OPTIX_ASSERT(optixDeviceContextSetLogCallback(optixDeviceContext,
                                                  contextOptions.logCallbackFunction,
                                                  contextOptions.logCallbackData,
                                                  contextOptions.logCallbackLevel));

#if (OPTIX_VERSION >= 90000)
    uint32_t clusterAccelFlags = OPTIX_DEVICE_PROPERTY_CLUSTER_ACCEL_FLAG_NONE;
    const OptixResult propertyResult =
        optixDeviceContextGetProperty(optixDeviceContext,
                                      OPTIX_DEVICE_PROPERTY_CLUSTER_ACCEL,
                                      &clusterAccelFlags,
                                      sizeof(clusterAccelFlags));
    if (propertyResult == OPTIX_SUCCESS)
    {
        supportsClusterBVHs =
            (clusterAccelFlags & OPTIX_DEVICE_PROPERTY_CLUSTER_ACCEL_FLAG_STANDARD) != 0;
    }
    else
#endif
    {
        supportsClusterBVHs = false;
    }

    if (!createRendererPipeline)
        return;
    CreatePipeline(MaterialBackend::OpenPbr);
}

bool CUDADevice::SupportsClusterBVHs() const
{
    return supportsClusterBVHs;
}

void CUDADevice::GetOptixCompileOptions(OptixModuleCompileOptions &moduleOptions,
                                        OptixPipelineCompileOptions &pipelineOptions)
{
    moduleOptions = {};
    moduleOptions.maxRegisterCount = OPTIX_COMPILE_DEFAULT_MAX_REGISTER_COUNT;
    moduleOptions.optLevel = OPTIX_COMPILE_OPTIMIZATION_DEFAULT;
    moduleOptions.debugLevel = OPTIX_COMPILE_DEBUG_LEVEL_MINIMAL;

    pipelineOptions = {};
    pipelineOptions.traversableGraphFlags = OPTIX_TRAVERSABLE_GRAPH_FLAG_ALLOW_ANY;
    pipelineOptions.usesMotionBlur = 0;
    pipelineOptions.usesPrimitiveTypeFlags = OPTIX_PRIMITIVE_TYPE_FLAGS_TRIANGLE;
    pipelineOptions.numPayloadValues = 2;
    pipelineOptions.numAttributeValues = 2;
    pipelineOptions.exceptionFlags = OPTIX_EXCEPTION_FLAG_NONE;
    pipelineOptions.pipelineLaunchParamsVariableName = "params";
    pipelineOptions.pipelineLaunchParamsSizeInBytes = sizeof(WavefrontParams);

    if (SupportsClusterBVHs())
    {
        pipelineOptions.allowClusteredGeometry = 1;
    }
}

void CUDADevice::CreatePipeline()
{
    CreatePipeline(MaterialBackend::OpenPbr);
}

void CUDADevice::DestroyPipeline()
{
    ContextScope contextScope(this);
    DestroyOptixPipeline();
}

bool CUDADevice::CreatePipeline(MaterialBackend backend)
{
    (void)backend;
    const char *ptxFilename = YBI_WAVEFRONT_KERNEL_PATH;
    std::ifstream input(ptxFilename, std::ios::binary);
    if (!input)
    {
        std::fprintf(stderr, "Failed to open PTX file: %s\n", ptxFilename);
        DestroyPipeline();
        return false;
    }
    const std::string ptx((std::istreambuf_iterator<char>(input)),
                          std::istreambuf_iterator<char>());
    if (ptx.empty())
    {
        std::fprintf(stderr, "PTX file is empty: %s\n", ptxFilename);
        DestroyPipeline();
        return false;
    }
    return CreateOptixPipeline(ptx);
}

bool CUDADevice::CreateOptixPipeline(const std::string &ptx)
{
    ContextScope contextScope(this);
    DestroyOptixPipeline();

    OptixModuleCompileOptions moduleCompileOptions = {};
    OptixPipelineCompileOptions pipelineCompileOptions = {};
    GetOptixCompileOptions(moduleCompileOptions, pipelineCompileOptions);

    char log[4096] = {};
    size_t logSize = sizeof(log);
    OptixResult result = CreateOptixModuleCompat(optixDeviceContext,
                                                 &moduleCompileOptions,
                                                 &pipelineCompileOptions,
                                                 ptx.c_str(),
                                                 ptx.size(),
                                                 log,
                                                 &logSize,
                                                 &optixPipeline.module);
    if (!CheckOptixResult(result, "optixModuleCreate", log, logSize))
    {
        DestroyOptixPipeline();
        return false;
    }

    OptixProgramGroupDesc programGroupDescs[PROGRAM_GROUP_COUNT] = {};
    OptixProgramGroupOptions programGroupOptions = {};
    static const char *entryNames[PROGRAM_GROUP_RAYGEN_END] = {
        "__raygen__intersect",
        "__raygen__shade_hit",
        "__raygen__direct_lighting",
        "__raygen__indirect_lighting",
    };
    for (int index = 0; index < PROGRAM_GROUP_RAYGEN_END; ++index)
    {
        OptixProgramGroupDesc &desc = programGroupDescs[index];
        desc.kind = OPTIX_PROGRAM_GROUP_KIND_RAYGEN;
        desc.raygen.module = optixPipeline.module;
        desc.raygen.entryFunctionName = entryNames[index];
    }

    OptixProgramGroupDesc &missDesc = programGroupDescs[PROGRAM_GROUP_MISS];
    missDesc.kind = OPTIX_PROGRAM_GROUP_KIND_MISS;
    missDesc.miss.module = optixPipeline.module;
    missDesc.miss.entryFunctionName = "__miss__kernel";

    OptixProgramGroupDesc &hitgroupDesc = programGroupDescs[PROGRAM_GROUP_CLOSEST_HIT];
    hitgroupDesc.kind = OPTIX_PROGRAM_GROUP_KIND_HITGROUP;
    hitgroupDesc.hitgroup.moduleCH = optixPipeline.module;
    hitgroupDesc.hitgroup.entryFunctionNameCH = "__closesthit__kernel";

    logSize = sizeof(log);
    result = optixProgramGroupCreate(optixDeviceContext,
                                     programGroupDescs,
                                     PROGRAM_GROUP_COUNT,
                                     &programGroupOptions,
                                     log,
                                     &logSize,
                                     optixPipeline.programGroups);
    if (!CheckOptixResult(result, "optixProgramGroupCreate", log, logSize))
    {
        DestroyOptixPipeline();
        return false;
    }

    for (int index = 0; index < PROGRAM_GROUP_COUNT; ++index)
    {
        result = optixSbtRecordPackHeader(optixPipeline.programGroups[index],
                                          &optixPipeline.sbtRecords[index]);
        if (!CheckOptixResult(result, "optixSbtRecordPackHeader"))
        {
            DestroyOptixPipeline();
            return false;
        }
    }

    const size_t recordsSize = sizeof(optixPipeline.sbtRecords);
    CUresult cudaResult = cuMemAlloc(&optixPipeline.sbtRecordsBuffer, recordsSize);
    if (!CheckCudaResult(cudaResult, "cuMemAlloc(SBT)"))
    {
        DestroyOptixPipeline();
        return false;
    }
    cudaResult =
        cuMemcpyHtoD(optixPipeline.sbtRecordsBuffer, optixPipeline.sbtRecords, recordsSize);
    if (!CheckCudaResult(cudaResult, "cuMemcpyHtoD(SBT)"))
    {
        DestroyOptixPipeline();
        return false;
    }

    if (!LinkOptixPipeline(pipelineCompileOptions))
    {
        DestroyOptixPipeline();
        return false;
    }
    return true;
}

bool CUDADevice::LinkOptixPipeline(const OptixPipelineCompileOptions &pipelineCompileOptions)
{
    DestroyLinkedOptixPipeline();

    std::vector<OptixProgramGroup> groups(optixPipeline.programGroups,
                                          optixPipeline.programGroups + PROGRAM_GROUP_COUNT);
    groups.insert(groups.end(),
                  optixPipeline.callableProgramGroups.begin(),
                  optixPipeline.callableProgramGroups.end());

    OptixPipelineLinkOptions pipelineLinkOptions = {};
    pipelineLinkOptions.maxTraceDepth = 1;
    char log[4096] = {};
    size_t logSize = sizeof(log);
    OptixResult result = optixPipelineCreate(optixDeviceContext,
                                             &pipelineCompileOptions,
                                             &pipelineLinkOptions,
                                             groups.data(),
                                             static_cast<unsigned int>(groups.size()),
                                             log,
                                             &logSize,
                                             &optixPipeline.pipeline);
    if (!CheckOptixResult(result, "optixPipelineCreate", log, logSize))
    {
        DestroyLinkedOptixPipeline();
        return false;
    }

    OptixStackSizes stackSizes = {};
    for (OptixProgramGroup group : groups)
    {
        result = AccumulateStackSizesCompat(group, &stackSizes, optixPipeline.pipeline);
        if (!CheckOptixResult(result, "optixUtilAccumulateStackSizes"))
        {
            DestroyLinkedOptixPipeline();
            return false;
        }
    }

    uint32_t directCallableStackSizeFromTraversal = 0;
    uint32_t directCallableStackSizeFromState = 0;
    uint32_t continuationStackSize = 0;
    const uint32_t maxDirectCallableDepth = optixPipeline.callableProgramGroups.empty() ? 0u : 1u;
    result = optixUtilComputeStackSizes(&stackSizes,
                                        1,
                                        0,
                                        maxDirectCallableDepth,
                                        &directCallableStackSizeFromTraversal,
                                        &directCallableStackSizeFromState,
                                        &continuationStackSize);
    if (!CheckOptixResult(result, "optixUtilComputeStackSizes"))
    {
        DestroyLinkedOptixPipeline();
        return false;
    }
    result = optixPipelineSetStackSize(optixPipeline.pipeline,
                                       directCallableStackSizeFromTraversal,
                                       directCallableStackSizeFromState,
                                       continuationStackSize,
                                       2);
    if (!CheckOptixResult(result, "optixPipelineSetStackSize"))
    {
        DestroyLinkedOptixPipeline();
        return false;
    }

    if (!optixPipeline.callableProgramGroups.empty())
    {
        std::vector<PipelineSbtRecord> records(optixPipeline.callableProgramGroups.size());
        for (size_t index = 0; index < records.size(); ++index)
        {
            result = optixSbtRecordPackHeader(optixPipeline.callableProgramGroups[index],
                                              &records[index]);
            if (!CheckOptixResult(result, "optixSbtRecordPackHeader(callable)"))
            {
                DestroyLinkedOptixPipeline();
                return false;
            }
        }

        const size_t recordsSize = sizeof(PipelineSbtRecord) * records.size();
        CUresult cudaResult = cuMemAlloc(&optixPipeline.callableRecords, recordsSize);
        if (!CheckCudaResult(cudaResult, "cuMemAlloc(callable SBT)"))
        {
            DestroyLinkedOptixPipeline();
            return false;
        }
        cudaResult = cuMemcpyHtoD(optixPipeline.callableRecords, records.data(), recordsSize);
        if (!CheckCudaResult(cudaResult, "cuMemcpyHtoD(callable SBT)"))
        {
            DestroyLinkedOptixPipeline();
            return false;
        }
    }
    return true;
}

void CUDADevice::DestroyLinkedOptixPipeline()
{
    if (optixPipeline.pipeline)
    {
        OptixResult result = optixPipelineDestroy(optixPipeline.pipeline);
        CheckOptixResult(result, "optixPipelineDestroy");
        optixPipeline.pipeline = nullptr;
    }

    if (optixPipeline.callableRecords)
    {
        CUresult result = cuMemFree(optixPipeline.callableRecords);
        CheckCudaResult(result, "cuMemFree(callable SBT)");
        optixPipeline.callableRecords = 0;
    }
}

void CUDADevice::DestroyOptixPipeline()
{
    DestroyLinkedOptixPipeline();

    if (optixPipeline.sbtRecordsBuffer)
    {
        CUresult result = cuMemFree(optixPipeline.sbtRecordsBuffer);
        CheckCudaResult(result, "cuMemFree(SBT)");
        optixPipeline.sbtRecordsBuffer = 0;
    }

    for (OptixProgramGroup programGroup : optixPipeline.callableProgramGroups)
    {
        OptixResult result = optixProgramGroupDestroy(programGroup);
        CheckOptixResult(result, "optixProgramGroupDestroy(callable)");
    }
    optixPipeline.callableProgramGroups.clear();

    for (int index = 0; index < PROGRAM_GROUP_COUNT; ++index)
    {
        if (optixPipeline.programGroups[index])
        {
            OptixResult result = optixProgramGroupDestroy(optixPipeline.programGroups[index]);
            CheckOptixResult(result, "optixProgramGroupDestroy");
            optixPipeline.programGroups[index] = nullptr;
        }
        optixPipeline.sbtRecords[index] = {};
    }

    if (optixPipeline.module)
    {
        OptixResult result = optixModuleDestroy(optixPipeline.module);
        CheckOptixResult(result, "optixModuleDestroy");
        optixPipeline.module = nullptr;
    }
    for (OptixModule module : optixPipeline.callableModules)
    {
        OptixResult result = optixModuleDestroy(module);
        CheckOptixResult(result, "optixModuleDestroy(callable)");
    }
    optixPipeline.callableModules.clear();
}

static OptixAccelBuildOptions GetDefaultBuildOptions(uint32_t numMotionKeys)
{
    OptixAccelBuildOptions options = {};
    options.buildFlags = OPTIX_BUILD_FLAG_ALLOW_COMPACTION | OPTIX_BUILD_FLAG_PREFER_FAST_TRACE;
    options.operation = OPTIX_BUILD_OPERATION_BUILD;
    options.motionOptions.numKeys = numMotionKeys;
    options.motionOptions.flags = 0;
    options.motionOptions.timeBegin = 0.f;
    options.motionOptions.timeEnd = 1.f;
    return options;
}

static void CopyTransformMatrix(const AffineTransform &transform, float out[12])
{
    memcpy(out, transform.m, sizeof(float) * 12);
}

static void SetOptixInstanceDefaults(OptixInstance &instance)
{
    memset(&instance, 0, sizeof(instance));
    instance.visibilityMask = OPTIX_VISIBILITY_ALL;
    instance.flags = OPTIX_INSTANCE_FLAG_NONE;
}

static constexpr unsigned int kInvalidInstanceId = ((1u << 24) - 1u);

static OptixTraversableHandle BuildOptixBVH(CUDADevice *cudaDevice,
                                            CUDAArena *deviceArena,
                                            OptixAccelBuildOptions buildOptions,
                                            const OptixBuildInput *buildInputs,
                                            uint32_t numBuildInputs)
{
    YBI_ASSERT(buildInputs);
    YBI_ASSERT(numBuildInputs > 0);

    OptixAccelBufferSizes sizes = {};
    OPTIX_ASSERT(optixAccelComputeMemoryUsage(
        cudaDevice->optixDeviceContext, &buildOptions, buildInputs, numBuildInputs, &sizes));

    cudaDevice->bvhTotalAllocated += sizes.outputSizeInBytes;

    const bool useCompaction = (buildOptions.buildFlags & OPTIX_BUILD_FLAG_ALLOW_COMPACTION) != 0;
    DeviceBuffer<std::byte> tempBuffer = deviceArena->PushArray<std::byte>(
        AlignUp(sizes.tempSizeInBytes + sizeof(uint64_t), sizeof(uint64_t)),
        OPTIX_ACCEL_BUFFER_BYTE_ALIGNMENT);

    DeviceBuffer<std::byte> outputBuffer =
        useCompaction ? deviceArena->PushArray<std::byte>(sizes.outputSizeInBytes,
                                                        OPTIX_ACCEL_BUFFER_BYTE_ALIGNMENT)
                      : cudaDevice->Alloc(sizes.outputSizeInBytes);

    uint64_t compactedSizeAddress =
        AlignUp(CUdeviceptr(tempBuffer.Data()) + sizes.tempSizeInBytes, sizeof(uint64_t));
    OptixAccelEmitDesc emittedProperties = {};
    emittedProperties.type = OPTIX_PROPERTY_TYPE_COMPACTED_SIZE;
    emittedProperties.result = CUdeviceptr(compactedSizeAddress);

    OptixTraversableHandle outputHandle = {};
    OPTIX_ASSERT(optixAccelBuild(cudaDevice->optixDeviceContext,
                                 0,
                                 &buildOptions,
                                 buildInputs,
                                 numBuildInputs,
                                 CUdeviceptr(tempBuffer.Data()),
                                 sizes.tempSizeInBytes,
                                 CUdeviceptr(outputBuffer.Data()),
                                 sizes.outputSizeInBytes,
                                 &outputHandle,
                                 useCompaction ? &emittedProperties : nullptr,
                                 useCompaction ? 1 : 0));

    CUDA_ASSERT(cuStreamSynchronize(0));

    if (useCompaction)
    {
        uint64_t compactedSize = sizes.outputSizeInBytes;
        CUDA_ASSERT(
            cuMemcpyDtoH(&compactedSize, CUdeviceptr(compactedSizeAddress), sizeof(uint64_t)));

        DeviceBuffer<std::byte> compactedOutputBuffer = cudaDevice->Alloc(compactedSize);

        OptixTraversableHandle tempOutputHandle = {};
        OPTIX_ASSERT(optixAccelCompact(cudaDevice->optixDeviceContext,
                                       0,
                                       outputHandle,
                                       CUdeviceptr(compactedOutputBuffer.Data()),
                                       compactedSize,
                                       &tempOutputHandle));
        CUDA_ASSERT(cuStreamSynchronize(0));
        outputHandle = tempOutputHandle;

        cudaDevice->bvhTotalAllocated -= sizes.outputSizeInBytes;
        cudaDevice->bvhTotalAllocated += compactedSize;

        // TEMP
        // cudaDevice->Free(compactedOutputBuffer);
    }

    return outputHandle;
}

static OptixBuildInput GetOptiXTriangleBuildInput(Arena *hostArena,
                                                  CUDAArena *deviceArena,
                                                  DataContainer *container,
                                                  Mesh &mesh,
                                                  uint32_t numMotionKeys)
{
    BufferHandle positionHandle = mesh.GetPositionsHandle();
    BufferHandle indexHandle = mesh.GetIndicesHandle();
    YBI_ASSERT(positionHandle.IsValid());
    YBI_ASSERT(indexHandle.IsValid());

    const Array<Vec3f> &positions = container->meshPositions[positionHandle.index];
    const Array<int> &indices = container->meshTriangleIndices[indexHandle.index];
    const uint32_t numVertices = static_cast<uint32_t>(positions.size());
    const uint32_t numIndices = static_cast<uint32_t>(indices.size());

    View<CUdeviceptr> vertexBuffers = hostArena->PushArray<CUdeviceptr>(numMotionKeys);
    View<Vec3f> hostVertices = hostArena->PushArray<Vec3f>(numVertices * numMotionKeys);

    DeviceBuffer<Vec3f> deviceVertices = deviceArena->PushArray<Vec3f>(numVertices * numMotionKeys);
    DeviceBuffer<int> deviceIndices = deviceArena->PushArray<int>(numIndices);

    for (uint32_t step = 0; step < numMotionKeys; step++)
    {
        CUdeviceptr dst = CUdeviceptr(deviceVertices.Data() + step * numVertices);
        vertexBuffers[step] = dst;

        memcpy(hostVertices.data() + step * numVertices,
               positions.data() + step * numVertices,
               sizeof(Vec3f) * numVertices);
    }

    CUDA_ASSERT(cuMemcpyHtoD(
        CUdeviceptr(deviceVertices.Data()), hostVertices.data(), deviceVertices.SizeInBytes()));
    CUDA_ASSERT(cuMemcpyHtoD(
        CUdeviceptr(deviceIndices.Data()), indices.data(), deviceIndices.SizeInBytes()));

    View<unsigned int> geometryFlags = hostArena->PushArray<unsigned int>(1);
    geometryFlags[0] = OPTIX_GEOMETRY_FLAG_DISABLE_ANYHIT;

    OptixBuildInput input = {};
    input.type = OPTIX_BUILD_INPUT_TYPE_TRIANGLES;
    OptixBuildInputTriangleArray &triangleArray = input.triangleArray;
    triangleArray = {};
    triangleArray.vertexBuffers = vertexBuffers.data();
    triangleArray.numVertices = numVertices;
    triangleArray.vertexFormat = OPTIX_VERTEX_FORMAT_FLOAT3;
    triangleArray.vertexStrideInBytes = sizeof(Vec3f);
    triangleArray.indexBuffer = CUdeviceptr(deviceIndices.Data());
    triangleArray.numIndexTriplets = numIndices / 3;
    triangleArray.indexFormat = OPTIX_INDICES_FORMAT_UNSIGNED_INT3;
    triangleArray.indexStrideInBytes = sizeof(int) * 3;
    triangleArray.preTransform = 0;
    triangleArray.transformFormat = OPTIX_TRANSFORM_FORMAT_NONE;
    triangleArray.flags = geometryFlags.data();
    triangleArray.numSbtRecords = 1;

    return input;
}

#if 0
static OptixBuildInput GetOptiXCurveBuildInput(CUDADevice *cudaDevice,
                                               HostMemoryArena &hostArena,
                                               CUDAMemoryArena &deviceArena,
                                               Curves &curves,
                                               uint32_t numMotionKeys,
                                               OptixAccelBuildOptions &options)
{
    (void)cudaDevice;
    (void)options;

    const uint32_t numVertices = (uint32_t)curves.GetNumVertices();
    const uint32_t totalNumSegments = (uint32_t)curves.GetNumSegments();
    const uint32_t totalNumVertices = numVertices * numMotionKeys;

    MemoryView<CUdeviceptr> vertexBuffers = hostArena.PushArray<CUdeviceptr>(numMotionKeys);
    MemoryView<CUdeviceptr> widthBuffers = hostArena.PushArray<CUdeviceptr>(numMotionKeys);
    MemoryView<uint32_t> indexBuffer = hostArena.PushArray<uint32_t>(totalNumSegments);
    MemoryView<Vec3> hostVertices = hostArena.PushArray<Vec3>(totalNumVertices);
    MemoryView<float> hostWidths = hostArena.PushArray<float>(totalNumVertices);

    DeviceMemoryView<Vec3> deviceVertices = deviceArena.PushArray<Vec3>(totalNumVertices);
    DeviceMemoryView<int> deviceIndices = deviceArena.PushArray<int>(totalNumSegments);
    DeviceMemoryView<float> deviceWidths = deviceArena.PushArray<float>(totalNumVertices);

    size_t segmentIndexOut = 0;
    for (size_t curveIndex = 0; curveIndex < curves.GetNumCurves(); curveIndex++)
    {
        const int segmentStart = curves.GetCurveKeyStart(curveIndex);
        const int numSegments = curves.GetCurveNumSegments(curveIndex);
        for (int segmentIndex = segmentStart; segmentIndex < segmentStart + numSegments;
             segmentIndex++, segmentIndexOut++)
        {
            indexBuffer[segmentIndexOut] = (uint32_t)segmentIndex;
        }
    }

    for (uint32_t step = 0; step < numMotionKeys; step++)
    {
        CUdeviceptr dst = (CUdeviceptr)(deviceVertices.data() + step * numVertices);
        vertexBuffers[step] = dst;

        dst = (CUdeviceptr)(deviceWidths.data() + step * numVertices);
        widthBuffers[step] = dst;

        const Array<Vec3> &positions = curves.GetVertices();
        const Array<float> &widths = curves.GetWidths();
        memcpy(hostVertices.data() + step * numVertices,
               positions.data(),
               sizeof(Vec3) * numVertices);
        memcpy(hostWidths.data() + step * numVertices, widths.data(), sizeof(float) * numVertices);
    }

    CUDA_ASSERT(cuMemcpyHtoD(
        CUdeviceptr(deviceVertices.data()), hostVertices.data(), deviceVertices.numBytes()));
    CUDA_ASSERT(cuMemcpyHtoD(
        CUdeviceptr(deviceIndices.data()), indexBuffer.data(), deviceIndices.numBytes()));
    CUDA_ASSERT(cuMemcpyHtoD(
        CUdeviceptr(deviceWidths.data()), hostWidths.data(), deviceWidths.numBytes()));

    OptixBuildInput input = {};
    input.type = OPTIX_BUILD_INPUT_TYPE_CURVES;
    OptixBuildInputCurveArray &curveArray = input.curveArray;
    curveArray = {};
    curveArray.curveType = OPTIX_PRIMITIVE_TYPE_ROUND_CUBIC_BSPLINE;
    curveArray.numPrimitives = totalNumSegments;
    curveArray.vertexBuffers = vertexBuffers.data();
    curveArray.numVertices = numVertices;
    curveArray.vertexStrideInBytes = sizeof(Vec3);
    curveArray.widthBuffers = widthBuffers.data();
    curveArray.widthStrideInBytes = sizeof(float);
    curveArray.normalBuffers = nullptr;
    curveArray.normalStrideInBytes = 0;
    curveArray.indexBuffer = CUdeviceptr(deviceIndices.data());
    curveArray.indexStrideInBytes = sizeof(int);
    curveArray.flag = OPTIX_GEOMETRY_FLAG_NONE;
    curveArray.primitiveIndexOffset = 0;
    curveArray.endcapFlags = OPTIX_CURVE_ENDCAP_DEFAULT;

    return input;
}
#endif

static OptixBuildInput GetOptiXInstanceBuildInput(CUDAArena *deviceArena,
                                                  View<OptixInstance> hostInstances)
{
    YBI_ASSERT(hostInstances.size() > 0);

    DeviceBuffer<OptixInstance> deviceInstances = deviceArena->PushArray<OptixInstance>(hostInstances.size());
    CUDA_ASSERT(cuMemcpyHtoD(
        CUdeviceptr(deviceInstances.Data()), hostInstances.data(), deviceInstances.SizeInBytes()));

    OptixBuildInput input = {};
    input.type = OPTIX_BUILD_INPUT_TYPE_INSTANCES;
    input.instanceArray.instances = CUdeviceptr(deviceInstances.Data());
    input.instanceArray.numInstances = (unsigned int)hostInstances.size();

    return input;
}

BVHHandle
CUDADevice::BuildTriangleBVH(DeviceArena *deviceArena, DataContainer *container, Mesh &mesh)
{
    ContextScope contextScope(this);
    ScratchArena scratch;

    DeviceArena::Scope scope(deviceArena);
    CUDAArena *cudaArena = reinterpret_cast<CUDAArena *>(deviceArena);

    const uint32_t numMotionKeys = 1;
    OptixAccelBuildOptions options = GetDefaultBuildOptions(numMotionKeys);

    OptixBuildInput buildInput =
        GetOptiXTriangleBuildInput(scratch.arena, cudaArena, container, mesh, numMotionKeys);
    OptixTraversableHandle handle = BuildOptixBVH(this, cudaArena, options, &buildInput, 1);

    return BVHHandle(handle);
}

BVHHandle CUDADevice::BuildInstancesBVH(DeviceArena *deviceArena,
                                        const Array<Instance> &instances,
                                        DataContainer *container,
                                        const Array<int> &objectToBVHIndex)
{
    YBI_ASSERT(instances.size());
    ContextScope contextScope(this);
    ScratchArena scratch;

    DeviceArena::Scope scope(deviceArena);
    CUDAArena *cudaArena = reinterpret_cast<CUDAArena *>(deviceArena);

    View<OptixInstance> sceneInstances = scratch.arena->PushArray<OptixInstance>(instances.size());

    std::vector<bool> lightMeshes(container->meshes.size(), false);
    for (const RectLight &light : container->rectLights)
    {
        YBI_ASSERT(light.meshIndex < lightMeshes.size());
        lightMeshes[light.meshIndex] = true;
    }

    for (size_t index = 0; index < instances.size(); index++)
    {
        const Instance &instance = instances[index];
        YBI_ASSERT(instance.handle.type != GeometryType::Scene);

        int bvhIndex = objectToBVHIndex[instance.handle.index];

        OptixInstance &optixInstance = sceneInstances[index];
        SetOptixInstanceDefaults(optixInstance);
        if (instance.handle.type == GeometryType::Mesh && lightMeshes[instance.handle.index])
        {
            optixInstance.visibilityMask = OPTIX_VISIBILITY_LIGHT;
        }
        CopyTransformMatrix(instance.transform, optixInstance.transform);
        optixInstance.instanceId = instance.handle.index;
        optixInstance.sbtOffset = 0;
        optixInstance.traversableHandle =
            (OptixTraversableHandle)container->GetBVHs(instance.handle.type)[bvhIndex]
                .handle.handle;
    }

    OptixBuildInput buildInput = GetOptiXInstanceBuildInput(cudaArena, sceneInstances);
    const uint32_t numMotionKeys = 1;
    OptixAccelBuildOptions options = GetDefaultBuildOptions(numMotionKeys);
    OptixTraversableHandle sceneHandle = BuildOptixBVH(this, cudaArena, options, &buildInput, 1);
    return BVHHandle(sceneHandle);
}

BVHHandle CUDADevice::BuildClusterBVH(DeviceArena *deviceArena,
                                      DataContainer *container,
                                      const Mesh &mesh,
                                      const MeshletClustersResult &clusterResult)
{
    (void)deviceArena;
#if (OPTIX_VERSION >= 90000)
    YBI_ASSERT(SupportsClusterBVHs());

    CUDA_ASSERT(cuCtxPushCurrent(cudaContext));
    YBI_ASSERT(clusterResult.meshlets.size() > 0);
    OptixTraversableHandle handle = {};
    BuildMeshCLAS(this, container, mesh, clusterResult, handle);

    CUDA_ASSERT(cuCtxPopCurrent(0));
    return BVHHandle(handle);
#else
    YBI_LOGFATAL("");
    return {};
#endif
}

size_t CUDADevice::GetBVHAllocationSize() const
{
    return bvhTotalAllocated;
}

void CUDADevice::StartProfiler() const
{
    cudaProfilerStart();
}

void CUDADevice::StopProfiler() const
{
    cudaProfilerStop();
}

#if 0
OptixTraversableHandle
BuildCurveGASFromCurves(CUDADevice *cudaDevice, HostMemoryArena &hostArena, Curves &curves)
{
    CUDA_ASSERT(cuCtxPushCurrent(cudaDevice->cudaContext));

    const uint32_t numMotionKeys = 1;
    OptixAccelBuildOptions options = GetDefaultBuildOptions(numMotionKeys);

    OptixBuildInput buildInput = GetOptiXCurveBuildInput(
        cudaDevice, hostArena, *cudaDevice->deviceArena, curves, numMotionKeys, options);
    OptixTraversableHandle handle =
        BuildOptixBVH(cudaDevice, *cudaDevice->deviceArena, options, &buildInput, 1);

    CUDA_ASSERT(cuCtxPopCurrent(0));
    return handle;
}

void CUDADevice::BuildBVH(Scene *scene)
{
    YBI_ASSERT(scene);

    CUDA_ASSERT(cuCtxPushCurrent(cudaContext));

    HostMemoryArena hostArena;
    std::vector<OptixInstance> sceneInstances;
    sceneInstances.reserve(scene->meshes.size() + scene->curves.size() + scene->instances.size());

    for (Mesh &mesh : scene->meshes)
    {
        if (mesh.positions.size() == 0 || mesh.indices.size() == 0)
        {
            continue;
        }

        OptixTraversableHandle meshHandle = BuildTriangleGASFromMesh(this, hostArena, mesh);
        if (!meshHandle)
        {
            continue;
        }

        OptixInstance optixInstance = {};
        SetOptixInstanceDefaults(optixInstance);
        CopyTransformMatrix(mesh.parentFromLocal, optixInstance.transform);
        YBI_ASSERT(mesh.refIndex != UINT32_MAX);
        YBI_ASSERT(mesh.refIndex < (1u << 24));
        optixInstance.instanceId = mesh.refIndex;
        optixInstance.sbtOffset = 0;
        optixInstance.traversableHandle = meshHandle;
        sceneInstances.push_back(optixInstance);

        hostArena.Clear();
        deviceArena->Clear();
    }

#if 0
    for (Curves &curves : scene->curves)
    {
        if (curves.GetNumVertices() == 0 || curves.GetNumCurves() == 0)
        {
            continue;
        }

        OptixTraversableHandle curvesHandle = BuildCurveGASFromCurves(this, hostArena, curves);
        if (!curvesHandle)
        {
            continue;
        }

        OptixInstance optixInstance = {};
        SetOptixInstanceDefaults(optixInstance);
        CopyTransformMatrix(curves.parentFromLocal, optixInstance.transform);
        optixInstance.instanceId = kInvalidInstanceId;
        optixInstance.sbtOffset = 0;
        optixInstance.traversableHandle = curvesHandle;
        sceneInstances.push_back(optixInstance);

        hostArena.Clear();
        deviceArena->Clear();
    }

    for (const Instance &instance : scene->instances)
    {
        YBI_ASSERT(instance.childSceneIndex < scene->childScenes.size());
        Scene *childScene = scene->childScenes[instance.childSceneIndex];
        YBI_ASSERT(childScene);
        YBI_ASSERT(childScene->bvhHandle != 0);

        OptixInstance optixInstance = {};
        SetOptixInstanceDefaults(optixInstance);
        CopyTransformMatrix(instance.parentFromLocal, optixInstance.transform);
        optixInstance.instanceId = kInvalidInstanceId;
        optixInstance.sbtOffset = 0;
        optixInstance.traversableHandle = (OptixTraversableHandle)childScene->bvhHandle;
        sceneInstances.push_back(optixInstance);
    }

    OptixTraversableHandle sceneHandle = {};
    if (!sceneInstances.empty())
    {
        MemoryView<OptixInstance> hostInstances =
            hostArena.PushArray<OptixInstance>(sceneInstances.size());
        memcpy(hostInstances.data(),
               sceneInstances.data(),
               sceneInstances.size() * sizeof(OptixInstance));
        OptixBuildInput buildInput = GetOptiXInstanceBuildInput(*deviceArena, hostInstances);
        const uint32_t numMotionKeys = 1;
        OptixAccelBuildOptions options = GetDefaultBuildOptions(numMotionKeys);
        sceneHandle = BuildOptixBVH(this, *deviceArena, options, &buildInput, 1);
    }
    scene->bvhHandle = (Scene::BVHHandle)sceneHandle;

    hostArena.Clear();
    deviceArena->Clear();
    CUDA_ASSERT(cuCtxPopCurrent(0));
}
#endif

#endif

#endif

} // namespace ybi
