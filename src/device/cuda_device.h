#pragma once

#include "cuda.h"
#include "cuda_runtime_api.h"
#include "optix_types.h"
#include "scene/scene.h"
#include "util/array.h"
#include "util/base.h"
#include "util/float3.h"

#include <cassert>
#include <cuda_runtime.h>
#include <optix_function_table_definition.h>
#include <optix_stubs.h>
#include <string>

YBI_NAMESPACE_BEGIN

#ifdef WITH_CUDA

#define CUDA_ASSERT(statement) \
    { \
        CUresult result = statement; \
        if (result != CUDA_SUCCESS) \
        { \
            const char *str, *name; \
            cuGetErrorString(result, &str); \
            cuGetErrorName(result, &name); \
            printf("CUDA Error (%s): %s in %s (%s:%d)\n", \
                   name, \
                   str, \
                   #statement, \
                   __FILE__, \
                   __LINE__); \
            assert(false); \
        } \
    }

#endif

#ifdef WITH_OPTIX

#define OPTIX_ASSERT(statement) \
    { \
        OptixResult result = statement; \
        if (result != OPTIX_SUCCESS) \
        { \
            const char *name = optixGetErrorString(result); \
            printf("Optix Error: %s in %s (%s:%d)\n", name, #statement, __FILE__, __LINE__); \
            assert(false); \
        } \
    }

struct ClusterAccelerationStructureLimits
{
    uint32_t maxTrianglesPerCluster;
    uint32_t maxVerticesPerCluster;
};

static OptixDeviceContext InitializeOptix(CUcontext cudaContext)
{
    OPTIX_ASSERT(optixInit());
    OptixDeviceContextOptions contextOptions = {};
    contextOptions.logCallbackFunction =
        [](unsigned int level, const char *tag, const char *message, void *cbdata) {
            std::string type = {};

            switch (level)
            {
                case 1:
                    type = "Fatal Error";
                    break;
                case 2:
                    type = "Error";
                    break;
                case 3:
                    type = "Warning";
                    break;
                case 4:
                    type = "Status";
                    break;
                default:
                    break;
            }

            // Print("Optix %S: %s\n", type, message);
        };
    contextOptions.logCallbackLevel = 4;
    contextOptions.validationMode = OPTIX_DEVICE_CONTEXT_VALIDATION_MODE_ALL;

    OptixDeviceContext optixDeviceContext;
    OPTIX_ASSERT(optixDeviceContextCreate(cudaContext, &contextOptions, &optixDeviceContext));
    OPTIX_ASSERT(optixDeviceContextSetLogCallback(optixDeviceContext,
                                                  contextOptions.logCallbackFunction,
                                                  contextOptions.logCallbackData,
                                                  contextOptions.logCallbackLevel));
    return optixDeviceContext;
}

// TODO: gpu kernels shouldn't be in the same file as host code
extern "C" __global__ void
Test(uint8_t *indexBuffer, float3 *vertexBuffer, ClusterAccelerationStructureLimits limits)
{
    uint32_t clusterId = 0;
    uint32_t triangleCount = 0;
    uint32_t vertexCount = 0;
    assert(triangleCount <= limits.maxTrianglesPerCluster);
    assert(vertexCount <= limits.maxVerticesPerCluster);

    OptixClusterAccelBuildInputTrianglesArgs args = {};
    args.clusterId = clusterId;
    args.clusterFlags = 0;
    args.triangleCount = triangleCount;
    args.vertexCount = vertexCount;
    args.positionTruncateBitCount = 0;
    args.indexFormat = OPTIX_CLUSTER_ACCEL_INDICES_FORMAT_8BIT;
    args.indexBuffer = (CUdeviceptr)indexBuffer;
    args.vertexBuffer = (CUdeviceptr)vertexBuffer;
}

struct CUDADevice
{
    CUcontext cudaContext;
    OptixDeviceContext optixDeviceContext;
    size_t totalAllocated;
    size_t bvhTotalAllocated;

    CUDADevice();
    void *Alloc(size_t size);
};

CUDADevice::CUDADevice() : totalAllocated(0), bvhTotalAllocated()
{
    cuInit(0);
    CUdevice device;
    cuDeviceGet(&device, 0);
    cuDevicePrimaryCtxRetain(&cudaContext, device);

    optixDeviceContext = InitializeOptix(cudaContext);
}

void *CUDADevice::Alloc(size_t size)
{
    assert(size != 0);
    totalAllocated += size;
    CUdeviceptr ptr;
    CUDA_ASSERT(cuMemAlloc(&ptr, size));
    return (void *)ptr;
}

static OptixTraversableHandle BuildOptixBVH(CUDADevice *cudaDevice,
                                            OptixAccelBuildOptions buildOptions,
                                            OptixBuildInput buildInput)
{
    OptixAccelBufferSizes sizes = {};
    OPTIX_ASSERT(optixAccelComputeMemoryUsage(
        cudaDevice->optixDeviceContext, &buildOptions, &buildInput, 1, &sizes));

    cudaDevice->bvhTotalAllocated += sizes.outputSizeInBytes;

    // Build BVH
    bool useCompaction = buildOptions.buildFlags & OPTIX_BUILD_FLAG_ALLOW_COMPACTION;
    void *tempBuffer =
        cudaDevice->Alloc(AlignUp(sizes.tempSizeInBytes + sizeof(uint64_t), sizeof(uint64_t)));
    void *outputBuffer = cudaDevice->Alloc(sizes.outputSizeInBytes);

    uint64_t compactedSizeAddress =
        AlignUp(CUdeviceptr(tempBuffer) + sizes.tempSizeInBytes, sizeof(uint64_t));
    OptixAccelEmitDesc emittedProperties = {};
    emittedProperties.type = OPTIX_PROPERTY_TYPE_COMPACTED_SIZE;
    emittedProperties.result = CUdeviceptr(compactedSizeAddress);
    OptixTraversableHandle outputHandle;
    OPTIX_ASSERT(optixAccelBuild(cudaDevice->optixDeviceContext,
                                 0,
                                 &buildOptions,
                                 &buildInput,
                                 1,
                                 CUdeviceptr(tempBuffer),
                                 sizes.tempSizeInBytes,
                                 CUdeviceptr(outputBuffer),
                                 sizes.outputSizeInBytes,
                                 &outputHandle,
                                 useCompaction ? &emittedProperties : 0,
                                 useCompaction ? 1 : 0));

    CUDA_ASSERT(cuStreamSynchronize(0));

    if (useCompaction)
    {
        uint64_t compactedSize = sizes.outputSizeInBytes;
        CUDA_ASSERT(
            cuMemcpyDtoH(&compactedSize, CUdeviceptr(compactedSizeAddress), sizeof(uint64_t)));

        if (compactedSize < sizes.outputSizeInBytes)
        {
            void *compactedOutputBuffer = cudaDevice->Alloc(compactedSize);

            OptixTraversableHandle tempOutputHandle;
            OPTIX_ASSERT(optixAccelCompact(cudaDevice->optixDeviceContext,
                                           0,
                                           outputHandle,
                                           CUdeviceptr(compactedOutputBuffer),
                                           compactedSize,
                                           &tempOutputHandle));
            CUDA_ASSERT(cuStreamSynchronize(0));
            outputHandle = tempOutputHandle;

            CUDA_ASSERT(cuMemFree(CUdeviceptr(outputBuffer)));
            cudaDevice->bvhTotalAllocated -= sizes.outputSizeInBytes;
            cudaDevice->bvhTotalAllocated += compactedSize;
        }
    }

    CUDA_ASSERT(cuMemFree(CUdeviceptr(tempBuffer)));
    return outputHandle;
}

// TODO: handle the case where cuda is enabled but optix isn't
void BuildBVH(CUDADevice *cudaDevice, Scene *scene)
{
    CUDA_ASSERT(cuCtxPushCurrent(cudaDevice->cudaContext));

    for (Mesh &mesh : scene->meshes)
    {
        const uint32_t numVertices = mesh.positions.size();
        const uint32_t numIndices = mesh.indices.size();
        uint32_t numMotionKeys = 1;

        OptixAccelBuildOptions options = {};
        options.buildFlags =
            OPTIX_BUILD_FLAG_ALLOW_COMPACTION | OPTIX_BUILD_FLAG_PREFER_FAST_TRACE;
        options.operation = OPTIX_BUILD_OPERATION_BUILD;
        options.motionOptions.numKeys = numMotionKeys;
        options.motionOptions.flags = 0;
        options.motionOptions.timeBegin = 0.f;
        options.motionOptions.timeEnd = 1.f;

        Array<CUdeviceptr> vertexBuffers(numMotionKeys);
        Array<float3> hostVertices(numVertices * numMotionKeys);

        size_t vertexSize = sizeof(float3) * numVertices * numMotionKeys;
        size_t indexSize = sizeof(int) * numIndices;
        float3 *deviceVertices = (float3 *)cudaDevice->Alloc(vertexSize);
        int *deviceIndices = (int *)cudaDevice->Alloc(indexSize);

        for (uint32_t step = 0; step < numMotionKeys; step++)
        {
            CUdeviceptr dst = (CUdeviceptr)(deviceVertices + step * numVertices);
            vertexBuffers[step] = dst;

            // TODO IMPORTANT: handle all motion blur data properly
            if (step > 0)
            {
                assert(0);
            }

            memcpy(hostVertices.data() + step * numVertices,
                   mesh.positions.data() + step * numVertices,
                   sizeof(float3) * numVertices);
        }
        CUDA_ASSERT(cuMemcpyHtoD(CUdeviceptr(deviceVertices), hostVertices.data(), vertexSize));
        CUDA_ASSERT(cuMemcpyHtoD(CUdeviceptr(deviceIndices), mesh.indices.data(), indexSize));

        unsigned int flags = OPTIX_GEOMETRY_FLAG_REQUIRE_SINGLE_ANYHIT_CALL;
        OptixBuildInput input = {};
        input.type = OPTIX_BUILD_INPUT_TYPE_TRIANGLES;
        OptixBuildInputTriangleArray &triangleArray = input.triangleArray;
        triangleArray = {};
        triangleArray.vertexBuffers = vertexBuffers.data();
        triangleArray.numVertices = numVertices;
        triangleArray.vertexFormat = OPTIX_VERTEX_FORMAT_FLOAT3;
        triangleArray.vertexStrideInBytes = 0;
        triangleArray.indexBuffer = CUdeviceptr(deviceIndices);
        triangleArray.numIndexTriplets = numIndices / 3;
        triangleArray.indexFormat = OPTIX_INDICES_FORMAT_UNSIGNED_INT3;
        triangleArray.indexStrideInBytes = 0;
        triangleArray.flags = &flags;
        triangleArray.numSbtRecords = 1;

        BuildOptixBVH(cudaDevice, options, input);
        CUDA_ASSERT(cuMemFree(CUdeviceptr(deviceVertices)));
        CUDA_ASSERT(cuMemFree(CUdeviceptr(deviceIndices)));
    }

    for (Curves &curve : scene->curves)
    {
        uint32_t numMotionKeys = 1;

        OptixAccelBuildOptions options = {};
        options.buildFlags =
            OPTIX_BUILD_FLAG_ALLOW_COMPACTION | OPTIX_BUILD_FLAG_PREFER_FAST_TRACE;
        options.operation = OPTIX_BUILD_OPERATION_BUILD;
        options.motionOptions.numKeys = numMotionKeys;
        options.motionOptions.flags = 0;
        options.motionOptions.timeBegin = 0.f;
        options.motionOptions.timeEnd = 1.f;

        Array<CUdeviceptr> vertexBuffers(numMotionKeys);
        Array<CUdeviceptr> widthBuffers(numMotionKeys);

        size_t totalNumSegments = curve.GetNumSegments();
        Array<uint32_t> indexBuffer(totalNumSegments);

        size_t bufferIndex = 0;
        for (size_t curveIndex = 0; curveIndex < curve.GetNumCurves(); curveIndex++)
        {
            int segmentStart = curve.GetCurveKeyStart(curveIndex);
            int numSegments = curve.GetCurveNumSegments(curveIndex);
            for (uint32_t segmentIndex = segmentStart; segmentIndex < segmentStart + numSegments;
                 segmentIndex++, bufferIndex++)
            {
                indexBuffer[bufferIndex] = segmentIndex;
            }
        }

        uint32_t numVertices = curve.GetNumVertices();
        uint32_t totalNumVertices = numVertices * numMotionKeys;
        Array<float3> hostVertices(totalNumVertices);
        Array<float> hostWidths(totalNumVertices);

        size_t vertexSize = sizeof(float3) * totalNumVertices;
        size_t indexSize = sizeof(int) * totalNumSegments;

        float3 *deviceVertices = (float3 *)cudaDevice->Alloc(vertexSize);
        int *deviceIndices = (int *)cudaDevice->Alloc(indexSize);
        float *deviceWidths = (float *)cudaDevice->Alloc(sizeof(float) * totalNumVertices);

        for (uint32_t step = 0; step < numMotionKeys; step++)
        {
            CUdeviceptr dst = (CUdeviceptr)(deviceVertices + step * numVertices);
            vertexBuffers[step] = dst;

            dst = (CUdeviceptr)(deviceWidths + step * numVertices);
            widthBuffers[step] = dst;

            const Array<float3> &positions = curve.GetVertices();
            const Array<float> &widths = curve.GetWidths();
            memcpy(hostVertices.data() + step * numVertices,
                   positions.data(),
                   sizeof(float3) * numVertices);
            memcpy(hostWidths.data() + step * numVertices,
                   widths.data(),
                   sizeof(float) * numVertices);
        }

        CUDA_ASSERT(cuMemcpyHtoD(CUdeviceptr(deviceVertices), hostVertices.data(), vertexSize));
        CUDA_ASSERT(cuMemcpyHtoD(CUdeviceptr(deviceIndices), indexBuffer.data(), indexSize));
        CUDA_ASSERT(cuMemcpyHtoD(
            CUdeviceptr(deviceWidths), hostWidths.data(), sizeof(float) * totalNumVertices));

        unsigned int flags = OPTIX_GEOMETRY_FLAG_REQUIRE_SINGLE_ANYHIT_CALL;

        OptixBuildInput input = {};
        input.type = OPTIX_BUILD_INPUT_TYPE_CURVES;
        OptixBuildInputCurveArray &curveArray = input.curveArray;
        curveArray = {};
        curveArray.curveType = OPTIX_PRIMITIVE_TYPE_ROUND_CUBIC_BSPLINE;
        curveArray.numPrimitives = totalNumSegments;
        curveArray.vertexBuffers = vertexBuffers.data();
        curveArray.numVertices = curve.GetNumVertices();
        curveArray.vertexStrideInBytes = 0;

        curveArray.widthBuffers = widthBuffers.data();
        curveArray.widthStrideInBytes = 0;
        curveArray.indexBuffer = CUdeviceptr(deviceIndices);
        curveArray.indexStrideInBytes = 0;
        curveArray.flag = flags;
        curveArray.primitiveIndexOffset = 0;
        curveArray.endcapFlags = OPTIX_CURVE_ENDCAP_DEFAULT;

        BuildOptixBVH(cudaDevice, options, input);
        CUDA_ASSERT(cuMemFree(CUdeviceptr(deviceVertices)));
        CUDA_ASSERT(cuMemFree(CUdeviceptr(deviceIndices)));
        CUDA_ASSERT(cuMemFree(CUdeviceptr(deviceWidths)));
    }

#if 0
#if (OPTIX_VERSION >= 90000)
    if (bvh->flags & BVHFlags::USE_CLUSTERS)
    {
        ClusterAccelerationStructureLimits limits;
        OPTIX_ASSERT(
            optixDeviceContextGetProperty(cudaDevice->optixDeviceContext,
                                          OPTIX_DEVICE_PROPERTY_LIMIT_MAX_CLUSTER_TRIANGLES,
                                          &limits.maxTrianglesPerCluster,
                                          sizeof(unsigned int)));
        OPTIX_ASSERT(
            optixDeviceContextGetProperty(cudaDevice->optixDeviceContext,
                                          OPTIX_DEVICE_PROPERTY_LIMIT_MAX_CLUSTER_VERTICES,
                                          &limits.maxVerticesPerCluster,
                                          sizeof(unsigned int)));
        printf("limits: %u %u\n", limits.maxVerticesPerCluster, limits.maxTrianglesPerCluster);

        OptixClusterAccelBuildInputTrianglesArgs *args;
        uint32_t *argsCount;
        uint32_t *outputSizes;
        void *tempBuffer;
        size_t tempBufferSize           = 0;
        uint32_t numClusters            = 0;
        uint32_t maxTrianglesPerCluster = 128;
        uint32_t maxVerticesPerCluster  = 256;

        OptixClusterAccelBuildModeDesc buildModeDesc = {};
        buildModeDesc.mode                          = OPTIX_CLUSTER_ACCEL_BUILD_MODE_GET_SIZES;
        buildModeDesc.getSize.outputSizesBuffer     = CUdeviceptr(outputSizes);
        buildModeDesc.getSize.tempBuffer            = CUdeviceptr(tempBuffer);
        buildModeDesc.getSize.tempBufferSizeInBytes = tempBufferSize;

        OptixClusterAccelBuildInput buildInput = {};
        buildInput.type            = OPTIX_CLUSTER_ACCEL_BUILD_TYPE_CLUSTERS_FROM_TRIANGLES;
        buildInput.triangles.flags = OPTIX_CLUSTER_ACCEL_BUILD_FLAG_PREFER_FAST_TRACE;
        buildInput.triangles.maxArgCount            = numClusters;
        buildInput.triangles.vertexFormat           = OPTIX_VERTEX_FORMAT_FLOAT3;
        buildInput.triangles.maxTriangleCountPerArg = maxTrianglesPerCluster;
        buildInput.triangles.maxVertexCountPerArg   = maxVerticesPerCluster;

        OPTIX_ASSERT(optixClusterAccelBuild(cudaDevice->optixDeviceContext,
                                            0,
                                            &buildModeDesc,
                                            &buildInput,
                                            CUdeviceptr(args),
                                            CUdeviceptr(argsCount),
                                            0));
    }
#else
    if (bvh->flags & BVHFlags::USE_CLUSTERS)
    {
        assert(0 && "Cluster Acceleration Structures are not supported. Please update OptiX.");
    }
#endif
#endif

    CUDA_ASSERT(cuCtxPopCurrent(0));
}

#endif

YBI_NAMESPACE_END
