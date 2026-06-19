#include "device/cuda_device.h"
#include "core/arena.h"
#include "core/assert.h"
#include "core/bit.h"
#include "core/thread_context.h"
#include "device/cuda_arena.h"
#include "device/cuda_device.h"
#include "scene/scene.h"
#include <cstddef>
#include <cstdio>
#include <cstring>
#include <optix_function_table_definition.h>
#include <optix_stack_size.h>
#include <optix_stubs.h>
#include <optix_types.h>
#include <string>
#include <vector>

namespace ybi
{

// Platform functions
#if defined(WITH_CUDA)

CUDADevice::ContextScope::ContextScope(CUDADevice *device)
{
    CUDA_ASSERT(cuCtxPushCurrent(device->cudaContext));
}

CUDADevice::ContextScope::~ContextScope()
{
    CUDA_ASSERT(cuCtxPopCurrent(0));
}

CUDADevice::CUDADevice() : bvhTotalAllocated(0)
{
    CUDA_ASSERT(cuInit(0));
    CUDA_ASSERT(cuDeviceGet(&device, 0));
    CUDA_ASSERT(cuDevicePrimaryCtxRetain(&cudaContext, device));

    InitializeOptix();
}

DeviceBuffer CUDADevice::Alloc(size_t numBytes)
{
    YBI_ASSERT(numBytes != 0);
    size_t size = numBytes;
    // totalAllocated += size;
    CUdeviceptr ptr;
    CUDA_ASSERT(cuMemAlloc(&ptr, size));

    return DeviceBuffer(uintptr_t(ptr), numBytes);
}

void CUDADevice::Free(DeviceBuffer &buffer)
{
    if (buffer.Data() == nullptr || buffer.SizeInBytes() == 0)
        return;
    // totalAllocated -= view.numBytes();
    CUDA_ASSERT(cuMemFree(CUdeviceptr(buffer.Data())));

    buffer = DeviceBuffer(std::nullptr_t());
}

std::unique_ptr<DeviceArena> CUDADevice::CreateArena()
{
    return std::make_unique<CUDAArena>();
}

#endif

// BVH
#if defined(WITH_CUDA) && defined(WITH_OPTIX)

#if 0
template <typename T>
struct alignas(OPTIX_SBT_RECORD_ALIGNMENT) PipelineSbtRecord
{
    char header[OPTIX_SBT_RECORD_HEADER_SIZE];
    T data;
};

struct EmptySbtData
{
};

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

#endif

void CUDADevice::InitializeOptix()
{
    const OptixResult initResult = optixInit();

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
}

#if 0

bool CUDADevice::CreateOptixPrimaryPipeline(const std::string &ptx)
{
    DestroyOptixPrimaryPipeline();

    OptixModuleCompileOptions moduleCompileOptions = {};
    moduleCompileOptions.maxRegisterCount = OPTIX_COMPILE_DEFAULT_MAX_REGISTER_COUNT;
    moduleCompileOptions.optLevel = OPTIX_COMPILE_OPTIMIZATION_DEFAULT;
    moduleCompileOptions.debugLevel = OPTIX_COMPILE_DEBUG_LEVEL_DEFAULT;

    OptixPipelineCompileOptions pipelineCompileOptions = {};
    pipelineCompileOptions.traversableGraphFlags = OPTIX_TRAVERSABLE_GRAPH_FLAG_ALLOW_ANY;
    pipelineCompileOptions.usesMotionBlur = 0;
    pipelineCompileOptions.usesPrimitiveTypeFlags = OPTIX_PRIMITIVE_TYPE_FLAGS_TRIANGLE |
                                                    OPTIX_PRIMITIVE_TYPE_FLAGS_ROUND_LINEAR |
                                                    OPTIX_PRIMITIVE_TYPE_FLAGS_ROUND_CUBIC_BSPLINE;
    pipelineCompileOptions.numPayloadValues = 2;
    pipelineCompileOptions.numAttributeValues = 4;
    pipelineCompileOptions.exceptionFlags = OPTIX_EXCEPTION_FLAG_NONE;
    pipelineCompileOptions.pipelineLaunchParamsVariableName = "params";

    char log[2048];
    size_t logSize = sizeof(log);
    if (!CheckOptixResult(CreateOptixModuleCompat(optixDeviceContext,
                                                  &moduleCompileOptions,
                                                  &pipelineCompileOptions,
                                                  ptx.c_str(),
                                                  ptx.size(),
                                                  log,
                                                  &logSize,
                                                  &optixPrimaryPipeline.module),
                          "CreateOptixModuleCompat"))
    {
        return false;
    }

    OptixBuiltinISOptions builtinISOptions = {};
    builtinISOptions.builtinISModuleType = OPTIX_PRIMITIVE_TYPE_ROUND_CUBIC_BSPLINE;
    builtinISOptions.usesMotionBlur = 0;
    builtinISOptions.buildFlags = OPTIX_BUILD_FLAG_PREFER_FAST_TRACE |
                                  OPTIX_BUILD_FLAG_ALLOW_COMPACTION |
                                  OPTIX_BUILD_FLAG_ALLOW_UPDATE;
    builtinISOptions.curveEndcapFlags = OPTIX_CURVE_ENDCAP_DEFAULT;
    if (!CheckOptixResult(optixBuiltinISModuleGet(optixDeviceContext,
                                                  &moduleCompileOptions,
                                                  &pipelineCompileOptions,
                                                  &builtinISOptions,
                                                  &optixPrimaryPipeline.curveModule),
                          "optixBuiltinISModuleGet"))
    {
        DestroyOptixPrimaryPipeline();
        return false;
    }

    OptixProgramGroupOptions programGroupOptions = {};
    OptixProgramGroupDesc raygenDesc = {};
    raygenDesc.kind = OPTIX_PROGRAM_GROUP_KIND_RAYGEN;
    raygenDesc.raygen.module = optixPrimaryPipeline.module;
    raygenDesc.raygen.entryFunctionName = "__raygen__primary";
    logSize = sizeof(log);
    if (!CheckOptixResult(optixProgramGroupCreate(optixDeviceContext,
                                                  &raygenDesc,
                                                  1,
                                                  &programGroupOptions,
                                                  log,
                                                  &logSize,
                                                  &optixPrimaryPipeline.raygenGroup),
                          "optixProgramGroupCreate(raygen)"))
    {
        DestroyOptixPrimaryPipeline();
        return false;
    }

    OptixProgramGroupDesc feedbackRaygenDesc = {};
    feedbackRaygenDesc.kind = OPTIX_PROGRAM_GROUP_KIND_RAYGEN;
    feedbackRaygenDesc.raygen.module = optixPrimaryPipeline.module;
    feedbackRaygenDesc.raygen.entryFunctionName = "__raygen__feedback";
    logSize = sizeof(log);
    if (!CheckOptixResult(optixProgramGroupCreate(optixDeviceContext,
                                                  &feedbackRaygenDesc,
                                                  1,
                                                  &programGroupOptions,
                                                  log,
                                                  &logSize,
                                                  &optixPrimaryPipeline.feedbackRaygenGroup),
                          "optixProgramGroupCreate(feedback raygen)"))
    {
        DestroyOptixPrimaryPipeline();
        return false;
    }

    OptixProgramGroupDesc missDesc = {};
    missDesc.kind = OPTIX_PROGRAM_GROUP_KIND_MISS;
    missDesc.miss.module = optixPrimaryPipeline.module;
    missDesc.miss.entryFunctionName = "__miss__primary";
    logSize = sizeof(log);
    if (!CheckOptixResult(optixProgramGroupCreate(optixDeviceContext,
                                                  &missDesc,
                                                  1,
                                                  &programGroupOptions,
                                                  log,
                                                  &logSize,
                                                  &optixPrimaryPipeline.missGroup),
                          "optixProgramGroupCreate(miss)"))
    {
        DestroyOptixPrimaryPipeline();
        return false;
    }

    OptixProgramGroupDesc hitgroupDesc = {};
    hitgroupDesc.kind = OPTIX_PROGRAM_GROUP_KIND_HITGROUP;
    hitgroupDesc.hitgroup.moduleCH = optixPrimaryPipeline.module;
    hitgroupDesc.hitgroup.entryFunctionNameCH = "__closesthit__primary";
    hitgroupDesc.hitgroup.moduleAH = optixPrimaryPipeline.module;
    hitgroupDesc.hitgroup.entryFunctionNameAH = "__anyhit__primary";
    hitgroupDesc.hitgroup.moduleIS = optixPrimaryPipeline.curveModule;
    hitgroupDesc.hitgroup.entryFunctionNameIS = nullptr;
    logSize = sizeof(log);
    if (!CheckOptixResult(optixProgramGroupCreate(optixDeviceContext,
                                                  &hitgroupDesc,
                                                  1,
                                                  &programGroupOptions,
                                                  log,
                                                  &logSize,
                                                  &optixPrimaryPipeline.hitgroupGroup),
                          "optixProgramGroupCreate(hitgroup)"))
    {
        DestroyOptixPrimaryPipeline();
        return false;
    }

    OptixProgramGroup groups[] = {
        optixPrimaryPipeline.raygenGroup,
        optixPrimaryPipeline.feedbackRaygenGroup,
        optixPrimaryPipeline.missGroup,
        optixPrimaryPipeline.hitgroupGroup,
    };

    OptixPipelineLinkOptions pipelineLinkOptions = {};
    pipelineLinkOptions.maxTraceDepth = 2;
    logSize = sizeof(log);
    if (!CheckOptixResult(optixPipelineCreate(optixDeviceContext,
                                              &pipelineCompileOptions,
                                              &pipelineLinkOptions,
                                              groups,
                                              4,
                                              log,
                                              &logSize,
                                              &optixPrimaryPipeline.pipeline),
                          "optixPipelineCreate"))
    {
        DestroyOptixPrimaryPipeline();
        return false;
    }

    OptixStackSizes stackSizes = {};
    if (!CheckOptixResult(AccumulateStackSizesCompat(optixPrimaryPipeline.raygenGroup,
                                                     &stackSizes,
                                                     optixPrimaryPipeline.pipeline),
                          "optixUtilAccumulateStackSizes(raygen)") ||
        !CheckOptixResult(AccumulateStackSizesCompat(optixPrimaryPipeline.feedbackRaygenGroup,
                                                     &stackSizes,
                                                     optixPrimaryPipeline.pipeline),
                          "optixUtilAccumulateStackSizes(feedback raygen)") ||
        !CheckOptixResult(AccumulateStackSizesCompat(optixPrimaryPipeline.missGroup,
                                                     &stackSizes,
                                                     optixPrimaryPipeline.pipeline),
                          "optixUtilAccumulateStackSizes(miss)") ||
        !CheckOptixResult(AccumulateStackSizesCompat(optixPrimaryPipeline.hitgroupGroup,
                                                     &stackSizes,
                                                     optixPrimaryPipeline.pipeline),
                          "optixUtilAccumulateStackSizes(hitgroup)"))
    {
        DestroyOptixPrimaryPipeline();
        return false;
    }

    uint32_t directCallableStackSizeFromTraversal = 0;
    uint32_t directCallableStackSizeFromState = 0;
    uint32_t continuationStackSize = 0;
    if (!CheckOptixResult(optixUtilComputeStackSizes(&stackSizes,
                                                     2,
                                                     0,
                                                     0,
                                                     &directCallableStackSizeFromTraversal,
                                                     &directCallableStackSizeFromState,
                                                     &continuationStackSize),
                          "optixUtilComputeStackSizes") ||
        !CheckOptixResult(optixPipelineSetStackSize(optixPrimaryPipeline.pipeline,
                                                    directCallableStackSizeFromTraversal,
                                                    directCallableStackSizeFromState,
                                                    continuationStackSize,
                                                    8),
                          "optixPipelineSetStackSize"))
    {
        DestroyOptixPrimaryPipeline();
        return false;
    }

    PipelineSbtRecord<EmptySbtData> raygenRecord = {};
    PipelineSbtRecord<EmptySbtData> feedbackRaygenRecord = {};
    PipelineSbtRecord<EmptySbtData> missRecord = {};
    PipelineSbtRecord<EmptySbtData> hitgroupRecord = {};
    if (!CheckOptixResult(
            optixSbtRecordPackHeader(optixPrimaryPipeline.raygenGroup, &raygenRecord),
            "optixSbtRecordPackHeader(raygen)") ||
        !CheckOptixResult(optixSbtRecordPackHeader(optixPrimaryPipeline.feedbackRaygenGroup,
                                                   &feedbackRaygenRecord),
                          "optixSbtRecordPackHeader(feedback raygen)") ||
        !CheckOptixResult(optixSbtRecordPackHeader(optixPrimaryPipeline.missGroup, &missRecord),
                          "optixSbtRecordPackHeader(miss)") ||
        !CheckOptixResult(
            optixSbtRecordPackHeader(optixPrimaryPipeline.hitgroupGroup, &hitgroupRecord),
            "optixSbtRecordPackHeader(hitgroup)"))
    {
        DestroyOptixPrimaryPipeline();
        return false;
    }

    CUDA_ASSERT(cuMemAlloc(&optixPrimaryPipeline.raygenRecordBuffer, sizeof(raygenRecord)));
    CUDA_ASSERT(cuMemAlloc(&optixPrimaryPipeline.feedbackRaygenRecordBuffer,
                           sizeof(feedbackRaygenRecord)));
    CUDA_ASSERT(cuMemAlloc(&optixPrimaryPipeline.missRecordBuffer, sizeof(missRecord)));
    CUDA_ASSERT(cuMemAlloc(&optixPrimaryPipeline.hitgroupRecordBuffer, sizeof(hitgroupRecord)));
    CUDA_ASSERT(cuMemcpyHtoD(
        optixPrimaryPipeline.raygenRecordBuffer, &raygenRecord, sizeof(raygenRecord)));
    CUDA_ASSERT(cuMemcpyHtoD(optixPrimaryPipeline.feedbackRaygenRecordBuffer,
                             &feedbackRaygenRecord,
                             sizeof(feedbackRaygenRecord)));
    CUDA_ASSERT(
        cuMemcpyHtoD(optixPrimaryPipeline.missRecordBuffer, &missRecord, sizeof(missRecord)));
    CUDA_ASSERT(cuMemcpyHtoD(
        optixPrimaryPipeline.hitgroupRecordBuffer, &hitgroupRecord, sizeof(hitgroupRecord)));

    optixPrimaryPipeline.sbt = {};
    optixPrimaryPipeline.sbt.raygenRecord = optixPrimaryPipeline.raygenRecordBuffer;
    optixPrimaryPipeline.sbt.missRecordBase = optixPrimaryPipeline.missRecordBuffer;
    optixPrimaryPipeline.sbt.missRecordStrideInBytes = sizeof(missRecord);
    optixPrimaryPipeline.sbt.missRecordCount = 1;
    optixPrimaryPipeline.sbt.hitgroupRecordBase = optixPrimaryPipeline.hitgroupRecordBuffer;
    optixPrimaryPipeline.sbt.hitgroupRecordStrideInBytes = sizeof(hitgroupRecord);
    optixPrimaryPipeline.sbt.hitgroupRecordCount = 1;
    return true;
}

void CUDADevice::DestroyOptixPrimaryPipeline()
{
    if (optixPrimaryPipeline.hitgroupRecordBuffer)
    {
        CUDA_ASSERT(cuMemFree(optixPrimaryPipeline.hitgroupRecordBuffer));
        optixPrimaryPipeline.hitgroupRecordBuffer = 0;
    }
    if (optixPrimaryPipeline.missRecordBuffer)
    {
        CUDA_ASSERT(cuMemFree(optixPrimaryPipeline.missRecordBuffer));
        optixPrimaryPipeline.missRecordBuffer = 0;
    }
    if (optixPrimaryPipeline.raygenRecordBuffer)
    {
        CUDA_ASSERT(cuMemFree(optixPrimaryPipeline.raygenRecordBuffer));
        optixPrimaryPipeline.raygenRecordBuffer = 0;
    }
    if (optixPrimaryPipeline.feedbackRaygenRecordBuffer)
    {
        CUDA_ASSERT(cuMemFree(optixPrimaryPipeline.feedbackRaygenRecordBuffer));
        optixPrimaryPipeline.feedbackRaygenRecordBuffer = 0;
    }
    if (optixPrimaryPipeline.pipeline)
    {
        OPTIX_ASSERT(optixPipelineDestroy(optixPrimaryPipeline.pipeline));
        optixPrimaryPipeline.pipeline = nullptr;
    }
    if (optixPrimaryPipeline.hitgroupGroup)
    {
        OPTIX_ASSERT(optixProgramGroupDestroy(optixPrimaryPipeline.hitgroupGroup));
        optixPrimaryPipeline.hitgroupGroup = nullptr;
    }
    if (optixPrimaryPipeline.missGroup)
    {
        OPTIX_ASSERT(optixProgramGroupDestroy(optixPrimaryPipeline.missGroup));
        optixPrimaryPipeline.missGroup = nullptr;
    }
    if (optixPrimaryPipeline.raygenGroup)
    {
        OPTIX_ASSERT(optixProgramGroupDestroy(optixPrimaryPipeline.raygenGroup));
        optixPrimaryPipeline.raygenGroup = nullptr;
    }
    if (optixPrimaryPipeline.feedbackRaygenGroup)
    {
        OPTIX_ASSERT(optixProgramGroupDestroy(optixPrimaryPipeline.feedbackRaygenGroup));
        optixPrimaryPipeline.feedbackRaygenGroup = nullptr;
    }
    if (optixPrimaryPipeline.curveModule)
    {
        OPTIX_ASSERT(optixModuleDestroy(optixPrimaryPipeline.curveModule));
        optixPrimaryPipeline.curveModule = nullptr;
    }
    if (optixPrimaryPipeline.module)
    {
        OPTIX_ASSERT(optixModuleDestroy(optixPrimaryPipeline.module));
        optixPrimaryPipeline.module = nullptr;
    }
    optixPrimaryPipeline.sbt = {};
}
#endif

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
    instance.visibilityMask = 0xFFu;
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
    DeviceBuffer tempBuffer = deviceArena->PushArray<uint8_t>(
        AlignUp(sizes.tempSizeInBytes + sizeof(uint64_t), sizeof(uint64_t)),
        OPTIX_ACCEL_BUFFER_BYTE_ALIGNMENT);

    DeviceBuffer outputBuffer =
        useCompaction ? deviceArena->PushArray<uint8_t>(sizes.outputSizeInBytes,
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

        DeviceBuffer compactedOutputBuffer = cudaDevice->Alloc(compactedSize);

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

    printf("%lli\n", cudaDevice->bvhTotalAllocated);

    return outputHandle;
}

static OptixBuildInput GetOptiXTriangleBuildInput(Arena *hostArena,
                                                  CUDAArena *deviceArena,
                                                  Mesh &mesh,
                                                  uint32_t numMotionKeys)
{
    const Array<Vec3f> &positions = mesh.GetPositions();
    const Array<int> &indices = mesh.GetIndices();
    const uint32_t numVertices = static_cast<uint32_t>(positions.size());
    const uint32_t numIndices = static_cast<uint32_t>(indices.size());

    View<CUdeviceptr> vertexBuffers = hostArena->PushArray<CUdeviceptr>(numMotionKeys);
    View<Vec3f> hostVertices = hostArena->PushArray<Vec3f>(numVertices * numMotionKeys);

    DeviceBuffer deviceVertices = deviceArena->PushArray<Vec3f>(numVertices * numMotionKeys);
    DeviceBuffer deviceIndices = deviceArena->PushArray<int>(numIndices);

    for (uint32_t step = 0; step < numMotionKeys; step++)
    {
        CUdeviceptr dst = CUdeviceptr(deviceVertices.As<Vec3f>().Get() + step * numVertices);
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

    DeviceBuffer deviceInstances = deviceArena->PushArray<OptixInstance>(hostInstances.size());
    CUDA_ASSERT(cuMemcpyHtoD(
        CUdeviceptr(deviceInstances.Data()), hostInstances.data(), deviceInstances.SizeInBytes()));

    OptixBuildInput input = {};
    input.type = OPTIX_BUILD_INPUT_TYPE_INSTANCES;
    input.instanceArray.instances = CUdeviceptr(deviceInstances.Data());
    input.instanceArray.numInstances = (unsigned int)hostInstances.size();

    return input;
}

BVHHandle CUDADevice::BuildTriangleBVH(DeviceArena *deviceArena, Mesh &mesh)
{
    ContextScope contextScope(this);
    ScratchArena scratch;

    DeviceArena::Scope scope(deviceArena);
    CUDAArena *cudaArena = reinterpret_cast<CUDAArena *>(deviceArena);

    const uint32_t numMotionKeys = 1;
    OptixAccelBuildOptions options = GetDefaultBuildOptions(numMotionKeys);

    OptixBuildInput buildInput =
        GetOptiXTriangleBuildInput(scratch.arena, cudaArena, mesh, numMotionKeys);
    OptixTraversableHandle handle = BuildOptixBVH(this, cudaArena, options, &buildInput, 1);

    return BVHHandle(handle);
}

BVHHandle CUDADevice::BuildInstancesBVH(DeviceArena *deviceArena,
                                        const Array<Instance> &instances,
                                        DataContainer *container)
{
    YBI_ASSERT(instances.size());
    ContextScope contextScope(this);
    ScratchArena scratch;

    DeviceArena::Scope scope(deviceArena);
    CUDAArena *cudaArena = reinterpret_cast<CUDAArena *>(deviceArena);

    View<OptixInstance> sceneInstances = scratch.arena->PushArray<OptixInstance>(instances.size());

    for (size_t index = 0; index < instances.size(); index++)
    {
        const Instance &instance = instances[index];
        YBI_ASSERT(instance.handle.type != GeometryType::Scene);

        OptixInstance &optixInstance = sceneInstances[index];
        SetOptixInstanceDefaults(optixInstance);
        CopyTransformMatrix(instance.transform, optixInstance.transform);
        optixInstance.instanceId = kInvalidInstanceId;
        optixInstance.sbtOffset = 0;
        optixInstance.traversableHandle =
            (OptixTraversableHandle)container->GetBVHs(instance.handle.type)[instance.handle.index]
                .handle.handle;
    }

    OptixBuildInput buildInput = GetOptiXInstanceBuildInput(cudaArena, sceneInstances);
    const uint32_t numMotionKeys = 1;
    OptixAccelBuildOptions options = GetDefaultBuildOptions(numMotionKeys);
    OptixTraversableHandle sceneHandle = BuildOptixBVH(this, cudaArena, options, &buildInput, 1);
    return BVHHandle(sceneHandle);
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
