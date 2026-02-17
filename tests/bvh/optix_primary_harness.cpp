#include "device/cuda_assert.h"
#include "io/usd/load.h"
#include "scene/scene.h"
#include "util/float3.h"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <limits>
#include <string>
#include <vector>

#include <cuda.h>
#include <optix.h>
#include <optix_function_table_definition.h>
#include <optix_stack_size.h>
#include <optix_stubs.h>

using namespace ybi;

#ifndef YBI_OPTIX_PRIMARY_PTX_PATH
#define YBI_OPTIX_PRIMARY_PTX_PATH ""
#endif

namespace
{
struct LaunchParams
{
    OptixTraversableHandle traversable;
    CUdeviceptr image;
    int width;
    int height;
    float3 cameraOrigin;
    float3 cameraU;
    float3 cameraV;
    float3 cameraW;
};

template <typename T>
struct alignas(OPTIX_SBT_RECORD_ALIGNMENT) SbtRecord
{
    char header[OPTIX_SBT_RECORD_HEADER_SIZE];
    T data;
};

struct EmptyData
{
};

using RaygenRecord = SbtRecord<EmptyData>;
using MissRecord = SbtRecord<EmptyData>;
using HitgroupRecord = SbtRecord<EmptyData>;

static float3 Cross(const float3 &a, const float3 &b)
{
    return make_float3(
        a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
}

static float3 Normalize(const float3 &v)
{
    const float lenSq = dot(v, v);
    if (lenSq <= 0.0f)
    {
        return make_float3(0.0f, 0.0f, 1.0f);
    }
    const float invLen = 1.0f / std::sqrt(lenSq);
    return v * invLen;
}

static bool SavePPM(const char *filePath, const std::vector<uint8_t> &rgba, int width, int height)
{
    std::ofstream out(filePath, std::ios::binary);
    if (!out.is_open())
    {
        return false;
    }

    out << "P6\n" << width << " " << height << "\n255\n";
    for (int y = 0; y < height; y++)
    {
        for (int x = 0; x < width; x++)
        {
            const size_t idx = static_cast<size_t>(y * width + x) * 4;
            out.put(static_cast<char>(rgba[idx + 0]));
            out.put(static_cast<char>(rgba[idx + 1]));
            out.put(static_cast<char>(rgba[idx + 2]));
        }
    }

    return out.good();
}

static std::string ReadTextFile(const std::string &path)
{
    std::ifstream input(path, std::ios::in | std::ios::binary);
    if (!input.is_open())
    {
        fprintf(stderr, "Failed to open PTX file: %s\n", path.c_str());
        std::abort();
    }

    return std::string((std::istreambuf_iterator<char>(input)), std::istreambuf_iterator<char>());
}

static void BuildFallbackMesh(Scene &scene)
{
    Mesh mesh;
    mesh.positions.Resize(8);
    mesh.indices.Resize(36);

    mesh.positions[0] = make_float3(-1.0f, -1.0f, -1.0f);
    mesh.positions[1] = make_float3(1.0f, -1.0f, -1.0f);
    mesh.positions[2] = make_float3(1.0f, 1.0f, -1.0f);
    mesh.positions[3] = make_float3(-1.0f, 1.0f, -1.0f);
    mesh.positions[4] = make_float3(-1.0f, -1.0f, 1.0f);
    mesh.positions[5] = make_float3(1.0f, -1.0f, 1.0f);
    mesh.positions[6] = make_float3(1.0f, 1.0f, 1.0f);
    mesh.positions[7] = make_float3(-1.0f, 1.0f, 1.0f);

    const int tris[] = {
        0, 1, 2, 0, 2, 3, // back
        4, 6, 5, 4, 7, 6, // front
        0, 4, 5, 0, 5, 1, // bottom
        3, 2, 6, 3, 6, 7, // top
        1, 5, 6, 1, 6, 2, // right
        0, 3, 7, 0, 7, 4  // left
    };

    for (int i = 0; i < 36; i++)
    {
        mesh.indices[i] = tris[i];
    }

    scene.meshes.emplace_back(std::move(mesh));
}
} // namespace

int main(int argc, char **argv)
{
    CUdevice device = 0;
    CUcontext cudaContext = nullptr;
    CUDA_ASSERT(cuInit(0));
    CUDA_ASSERT(cuDeviceGet(&device, 0));
    CUDA_ASSERT(cuDevicePrimaryCtxRetain(&cudaContext, device));
    CUDA_ASSERT(cuCtxSetCurrent(cudaContext));

    OPTIX_ASSERT(optixInit());

    OptixDeviceContextOptions contextOptions = {};
    contextOptions.logCallbackFunction = [](unsigned int level,
                                            const char *tag,
                                            const char *message,
                                            void *cbdata) {
        (void)cbdata;
        printf("[OptiX][%u][%s] %s\n", level, tag ? tag : "", message ? message : "");
    };
    contextOptions.logCallbackLevel = 4;
#ifdef YBI_DEBUG
    contextOptions.validationMode = OPTIX_DEVICE_CONTEXT_VALIDATION_MODE_ALL;
#endif

    OptixDeviceContext optixContext = nullptr;
    OPTIX_ASSERT(optixDeviceContextCreate(cudaContext, &contextOptions, &optixContext));

    bool loadUsd = false;
    for (int argIndex = 1; argIndex < argc; argIndex++)
    {
        if (std::string(argv[argIndex]) == "--usd")
        {
            loadUsd = true;
        }
    }

    Scene scene;
    if (loadUsd)
    {
        LoadUSDScene(&scene);
    }
    if (scene.meshes.empty())
    {
        printf("Using fallback cube mesh.\n");
        BuildFallbackMesh(scene);
    }

    Mesh &mesh = scene.meshes[0];
    const uint32_t numVertices = static_cast<uint32_t>(mesh.positions.size());
    const uint32_t numIndices = static_cast<uint32_t>(mesh.indices.size());
    if (numVertices == 0 || numIndices == 0)
    {
        fprintf(stderr, "Mesh has no geometry.\n");
        return 1;
    }

    float3 boundsMin = make_float3(
        std::numeric_limits<float>::max(),
        std::numeric_limits<float>::max(),
        std::numeric_limits<float>::max());
    float3 boundsMax = make_float3(
        -std::numeric_limits<float>::max(),
        -std::numeric_limits<float>::max(),
        -std::numeric_limits<float>::max());
    for (uint32_t i = 0; i < numVertices; i++)
    {
        const float3 p = mesh.positions[i];
        boundsMin.x = std::min(boundsMin.x, p.x);
        boundsMin.y = std::min(boundsMin.y, p.y);
        boundsMin.z = std::min(boundsMin.z, p.z);
        boundsMax.x = std::max(boundsMax.x, p.x);
        boundsMax.y = std::max(boundsMax.y, p.y);
        boundsMax.z = std::max(boundsMax.z, p.z);
    }

    const float3 boundsCenter = (boundsMin + boundsMax) * 0.5f;
    const float3 boundsExtent = boundsMax - boundsMin;
    const float diagonal = std::max(0.001f, length(boundsExtent));

    CUdeviceptr vertexBuffer = 0;
    CUdeviceptr indexBuffer = 0;
    CUDA_ASSERT(cuMemAlloc(&vertexBuffer, sizeof(float3) * numVertices));
    CUDA_ASSERT(cuMemAlloc(&indexBuffer, sizeof(int) * numIndices));
    CUDA_ASSERT(cuMemcpyHtoD(vertexBuffer, mesh.positions.data(), sizeof(float3) * numVertices));
    CUDA_ASSERT(cuMemcpyHtoD(indexBuffer, mesh.indices.data(), sizeof(int) * numIndices));

    uint32_t triangleInputFlags[1] = {OPTIX_GEOMETRY_FLAG_DISABLE_ANYHIT};
    OptixBuildInput buildInput = {};
    buildInput.type = OPTIX_BUILD_INPUT_TYPE_TRIANGLES;
    buildInput.triangleArray.vertexFormat = OPTIX_VERTEX_FORMAT_FLOAT3;
    buildInput.triangleArray.vertexStrideInBytes = sizeof(float3);
    buildInput.triangleArray.numVertices = numVertices;
    buildInput.triangleArray.vertexBuffers = &vertexBuffer;
    buildInput.triangleArray.indexFormat = OPTIX_INDICES_FORMAT_UNSIGNED_INT3;
    buildInput.triangleArray.indexStrideInBytes = sizeof(int) * 3;
    buildInput.triangleArray.numIndexTriplets = numIndices / 3;
    buildInput.triangleArray.indexBuffer = indexBuffer;
    buildInput.triangleArray.flags = triangleInputFlags;
    buildInput.triangleArray.numSbtRecords = 1;

    OptixAccelBuildOptions accelOptions = {};
    accelOptions.buildFlags = OPTIX_BUILD_FLAG_ALLOW_COMPACTION | OPTIX_BUILD_FLAG_PREFER_FAST_TRACE;
    accelOptions.operation = OPTIX_BUILD_OPERATION_BUILD;

    OptixAccelBufferSizes gasBufferSizes = {};
    OPTIX_ASSERT(optixAccelComputeMemoryUsage(
        optixContext, &accelOptions, &buildInput, 1, &gasBufferSizes));

    CUdeviceptr tempBuffer = 0;
    CUdeviceptr outputBuffer = 0;
    CUDA_ASSERT(cuMemAlloc(&tempBuffer, gasBufferSizes.tempSizeInBytes));
    CUDA_ASSERT(cuMemAlloc(&outputBuffer, gasBufferSizes.outputSizeInBytes));

    CUdeviceptr compactedSizeBuffer = 0;
    CUDA_ASSERT(cuMemAlloc(&compactedSizeBuffer, sizeof(uint64_t)));
    OptixAccelEmitDesc emitProperty = {};
    emitProperty.type = OPTIX_PROPERTY_TYPE_COMPACTED_SIZE;
    emitProperty.result = compactedSizeBuffer;

    OptixTraversableHandle gasHandle = 0;
    OPTIX_ASSERT(optixAccelBuild(optixContext,
                                 0,
                                 &accelOptions,
                                 &buildInput,
                                 1,
                                 tempBuffer,
                                 gasBufferSizes.tempSizeInBytes,
                                 outputBuffer,
                                 gasBufferSizes.outputSizeInBytes,
                                 &gasHandle,
                                 &emitProperty,
                                 1));
    CUDA_ASSERT(cuStreamSynchronize(0));

    uint64_t compactedSize = 0;
    CUDA_ASSERT(cuMemcpyDtoH(&compactedSize, compactedSizeBuffer, sizeof(uint64_t)));
    CUdeviceptr gasBuffer = outputBuffer;
    if (compactedSize > 0 && compactedSize < gasBufferSizes.outputSizeInBytes)
    {
        CUdeviceptr compactedBuffer = 0;
        CUDA_ASSERT(cuMemAlloc(&compactedBuffer, compactedSize));
        OPTIX_ASSERT(optixAccelCompact(
            optixContext, 0, gasHandle, compactedBuffer, compactedSize, &gasHandle));
        CUDA_ASSERT(cuStreamSynchronize(0));
        CUDA_ASSERT(cuMemFree(outputBuffer));
        gasBuffer = compactedBuffer;
    }

    const std::string ptx = ReadTextFile(YBI_OPTIX_PRIMARY_PTX_PATH);

    OptixModuleCompileOptions moduleCompileOptions = {};
    moduleCompileOptions.maxRegisterCount = OPTIX_COMPILE_DEFAULT_MAX_REGISTER_COUNT;
    moduleCompileOptions.optLevel = OPTIX_COMPILE_OPTIMIZATION_DEFAULT;
    moduleCompileOptions.debugLevel = OPTIX_COMPILE_DEBUG_LEVEL_DEFAULT;

    OptixPipelineCompileOptions pipelineCompileOptions = {};
    pipelineCompileOptions.traversableGraphFlags = OPTIX_TRAVERSABLE_GRAPH_FLAG_ALLOW_SINGLE_GAS;
    pipelineCompileOptions.usesMotionBlur = 0;
    pipelineCompileOptions.numPayloadValues = 1;
    pipelineCompileOptions.numAttributeValues = 2;
    pipelineCompileOptions.exceptionFlags = OPTIX_EXCEPTION_FLAG_NONE;
    pipelineCompileOptions.pipelineLaunchParamsVariableName = "params";

    char log[2048];
    size_t logSize = sizeof(log);
    OptixModule module = nullptr;
    OPTIX_ASSERT(optixModuleCreate(optixContext,
                                   &moduleCompileOptions,
                                   &pipelineCompileOptions,
                                   ptx.c_str(),
                                   ptx.size(),
                                   log,
                                   &logSize,
                                   &module));
    if (logSize > 1)
    {
        printf("%s\n", log);
    }

    OptixProgramGroupOptions programGroupOptions = {};

    OptixProgramGroupDesc raygenDesc = {};
    raygenDesc.kind = OPTIX_PROGRAM_GROUP_KIND_RAYGEN;
    raygenDesc.raygen.module = module;
    raygenDesc.raygen.entryFunctionName = "__raygen__primary";
    OptixProgramGroup raygenGroup = nullptr;
    logSize = sizeof(log);
    OPTIX_ASSERT(
        optixProgramGroupCreate(
            optixContext, &raygenDesc, 1, &programGroupOptions, log, &logSize, &raygenGroup));

    OptixProgramGroupDesc missDesc = {};
    missDesc.kind = OPTIX_PROGRAM_GROUP_KIND_MISS;
    missDesc.miss.module = module;
    missDesc.miss.entryFunctionName = "__miss__primary";
    OptixProgramGroup missGroup = nullptr;
    logSize = sizeof(log);
    OPTIX_ASSERT(
        optixProgramGroupCreate(
            optixContext, &missDesc, 1, &programGroupOptions, log, &logSize, &missGroup));

    OptixProgramGroupDesc hitgroupDesc = {};
    hitgroupDesc.kind = OPTIX_PROGRAM_GROUP_KIND_HITGROUP;
    hitgroupDesc.hitgroup.moduleCH = module;
    hitgroupDesc.hitgroup.entryFunctionNameCH = "__closesthit__primary";
    OptixProgramGroup hitgroupGroup = nullptr;
    logSize = sizeof(log);
    OPTIX_ASSERT(optixProgramGroupCreate(optixContext,
                                         &hitgroupDesc,
                                         1,
                                         &programGroupOptions,
                                         log,
                                         &logSize,
                                         &hitgroupGroup));

    OptixProgramGroup groups[] = {raygenGroup, missGroup, hitgroupGroup};
    OptixPipelineLinkOptions pipelineLinkOptions = {};
    pipelineLinkOptions.maxTraceDepth = 1;

    OptixPipeline pipeline = nullptr;
    logSize = sizeof(log);
    OPTIX_ASSERT(optixPipelineCreate(optixContext,
                                     &pipelineCompileOptions,
                                     &pipelineLinkOptions,
                                     groups,
                                     3,
                                     log,
                                     &logSize,
                                     &pipeline));
    if (logSize > 1)
    {
        printf("%s\n", log);
    }

    OptixStackSizes stackSizes = {};
    OPTIX_ASSERT(optixUtilAccumulateStackSizes(raygenGroup, &stackSizes, pipeline));
    OPTIX_ASSERT(optixUtilAccumulateStackSizes(missGroup, &stackSizes, pipeline));
    OPTIX_ASSERT(optixUtilAccumulateStackSizes(hitgroupGroup, &stackSizes, pipeline));
    uint32_t directCallableStackSizeFromTraversal = 0;
    uint32_t directCallableStackSizeFromState = 0;
    uint32_t continuationStackSize = 0;
    OPTIX_ASSERT(optixUtilComputeStackSizes(&stackSizes,
                                            1,
                                            0,
                                            0,
                                            &directCallableStackSizeFromTraversal,
                                            &directCallableStackSizeFromState,
                                            &continuationStackSize));
    OPTIX_ASSERT(optixPipelineSetStackSize(pipeline,
                                           directCallableStackSizeFromTraversal,
                                           directCallableStackSizeFromState,
                                           continuationStackSize,
                                           1));

    RaygenRecord raygenRecord = {};
    MissRecord missRecord = {};
    HitgroupRecord hitgroupRecord = {};
    OPTIX_ASSERT(optixSbtRecordPackHeader(raygenGroup, &raygenRecord));
    OPTIX_ASSERT(optixSbtRecordPackHeader(missGroup, &missRecord));
    OPTIX_ASSERT(optixSbtRecordPackHeader(hitgroupGroup, &hitgroupRecord));

    CUdeviceptr raygenRecordBuffer = 0;
    CUdeviceptr missRecordBuffer = 0;
    CUdeviceptr hitgroupRecordBuffer = 0;
    CUDA_ASSERT(cuMemAlloc(&raygenRecordBuffer, sizeof(RaygenRecord)));
    CUDA_ASSERT(cuMemAlloc(&missRecordBuffer, sizeof(MissRecord)));
    CUDA_ASSERT(cuMemAlloc(&hitgroupRecordBuffer, sizeof(HitgroupRecord)));
    CUDA_ASSERT(cuMemcpyHtoD(raygenRecordBuffer, &raygenRecord, sizeof(RaygenRecord)));
    CUDA_ASSERT(cuMemcpyHtoD(missRecordBuffer, &missRecord, sizeof(MissRecord)));
    CUDA_ASSERT(cuMemcpyHtoD(hitgroupRecordBuffer, &hitgroupRecord, sizeof(HitgroupRecord)));

    OptixShaderBindingTable sbt = {};
    sbt.raygenRecord = raygenRecordBuffer;
    sbt.missRecordBase = missRecordBuffer;
    sbt.missRecordStrideInBytes = sizeof(MissRecord);
    sbt.missRecordCount = 1;
    sbt.hitgroupRecordBase = hitgroupRecordBuffer;
    sbt.hitgroupRecordStrideInBytes = sizeof(HitgroupRecord);
    sbt.hitgroupRecordCount = 1;

    const int width = 1280;
    const int height = 720;
    const size_t imageSize = static_cast<size_t>(width) * static_cast<size_t>(height) * 4;
    CUdeviceptr imageBuffer = 0;
    CUDA_ASSERT(cuMemAlloc(&imageBuffer, imageSize));

    const float3 eye = boundsCenter + make_float3(0.0f, 0.0f, 2.5f * diagonal);
    const float3 target = boundsCenter;
    const float3 worldUp = make_float3(0.0f, 1.0f, 0.0f);
    const float3 forward = Normalize(target - eye);
    const float3 right = Normalize(Cross(forward, worldUp));
    const float3 up = Normalize(Cross(right, forward));
    const float aspect = static_cast<float>(width) / static_cast<float>(height);
    const float fovY = 45.0f * 3.14159265358979323846f / 180.0f;
    const float tanHalfFov = std::tan(fovY * 0.5f);

    LaunchParams params = {};
    params.traversable = gasHandle;
    params.image = imageBuffer;
    params.width = width;
    params.height = height;
    params.cameraOrigin = eye;
    params.cameraU = right * (aspect * tanHalfFov);
    params.cameraV = up * tanHalfFov;
    params.cameraW = forward;

    CUdeviceptr paramsBuffer = 0;
    CUDA_ASSERT(cuMemAlloc(&paramsBuffer, sizeof(LaunchParams)));
    CUDA_ASSERT(cuMemcpyHtoD(paramsBuffer, &params, sizeof(LaunchParams)));

    OPTIX_ASSERT(optixLaunch(
        pipeline, 0, paramsBuffer, sizeof(LaunchParams), &sbt, width, height, 1));
    CUDA_ASSERT(cuStreamSynchronize(0));

    std::vector<uint8_t> hostImage(imageSize);
    CUDA_ASSERT(cuMemcpyDtoH(hostImage.data(), imageBuffer, imageSize));

    const char *outputFile = "optix_primary.ppm";
    if (!SavePPM(outputFile, hostImage, width, height))
    {
        fprintf(stderr, "Failed to write image: %s\n", outputFile);
        return 1;
    }

    printf("Wrote %s (%dx%d)\n", outputFile, width, height);
    printf("Triangles: %u, GAS bytes: %llu\n",
           numIndices / 3,
           static_cast<unsigned long long>(compactedSize > 0 ? compactedSize : gasBufferSizes.outputSizeInBytes));

    CUDA_ASSERT(cuMemFree(paramsBuffer));
    CUDA_ASSERT(cuMemFree(imageBuffer));
    CUDA_ASSERT(cuMemFree(hitgroupRecordBuffer));
    CUDA_ASSERT(cuMemFree(missRecordBuffer));
    CUDA_ASSERT(cuMemFree(raygenRecordBuffer));
    CUDA_ASSERT(cuMemFree(tempBuffer));
    CUDA_ASSERT(cuMemFree(compactedSizeBuffer));
    CUDA_ASSERT(cuMemFree(gasBuffer));
    CUDA_ASSERT(cuMemFree(indexBuffer));
    CUDA_ASSERT(cuMemFree(vertexBuffer));

    OPTIX_ASSERT(optixPipelineDestroy(pipeline));
    OPTIX_ASSERT(optixProgramGroupDestroy(hitgroupGroup));
    OPTIX_ASSERT(optixProgramGroupDestroy(missGroup));
    OPTIX_ASSERT(optixProgramGroupDestroy(raygenGroup));
    OPTIX_ASSERT(optixModuleDestroy(module));
    OPTIX_ASSERT(optixDeviceContextDestroy(optixContext));

    CUDA_ASSERT(cuDevicePrimaryCtxRelease(device));
    return 0;
}
