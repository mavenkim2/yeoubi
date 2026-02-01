#pragma once

#include "optix_types.h"
#include <cassert>
#include <cuda_runtime.h>
#include <optix_function_table_definition.h>
#include <optix_stubs.h>
#include <string>

YBI_NAMESPACE_BEGIN

#ifdef WITH_OPTIX

struct ClusterAccelerationStructureLimits
{
    uint32_t maxTrianglesPerCluster;
    uint32_t maxVerticesPerCluster;
};

static OptixDeviceContext InitializeOptix(CUcontext cudaContext)
{
    optixInit();
    OptixDeviceContextOptions contextOptions = {};
    contextOptions.logCallbackFunction       = [](unsigned int level, const char *tag,
                                            const char *message, void *cbdata) {
        std::string type = {};

        switch (level)
        {
            case 1: type = "Fatal Error"; break;
            case 2: type = "Error"; break;
            case 3: type = "Warning"; break;
            case 4: type = "Status"; break;
            default: break;
        }

        // Print("Optix %S: %s\n", type, message);
    };
    contextOptions.logCallbackLevel = 4;
    contextOptions.validationMode   = OPTIX_DEVICE_CONTEXT_VALIDATION_MODE_ALL;

    OptixDeviceContext optixDeviceContext;
    optixDeviceContextCreate(cudaContext, &contextOptions, &optixDeviceContext);
    optixDeviceContextSetLogCallback(optixDeviceContext, contextOptions.logCallbackFunction,
                                     contextOptions.logCallbackData,
                                     contextOptions.logCallbackLevel);
    return optixDeviceContext;
}

extern "C" __global__ void Test(uint8_t *indexBuffer, float3 *vertexBuffer,
                                ClusterAccelerationStructureLimits limits)
{
    uint32_t clusterId     = 0;
    uint32_t triangleCount = 0;
    uint32_t vertexCount   = 0;
    assert(triangleCount <= limits.maxTrianglesPerCluster);
    assert(vertexCount <= limits.maxVerticesPerCluster);

    OptixClusterAccelBuildInputTrianglesArgs args = {};
    args.clusterId                                = clusterId;
    args.clusterFlags                             = 0;
    args.triangleCount                            = triangleCount;
    args.vertexCount                              = vertexCount;
    args.positionTruncateBitCount                 = 0;
    args.indexFormat                              = OPTIX_CLUSTER_ACCEL_INDICES_FORMAT_8BIT;
    args.indexBuffer                              = (CUdeviceptr)indexBuffer;
    args.vertexBuffer                             = (CUdeviceptr)vertexBuffer;
}

static void BuildBVH()
{
    enum BVHFlags : int
    {
        USE_CLUSTERS = (1u << 0u),
    };
    struct BVH
    {
        BVHFlags flags;
    };

    cuInit(0);
    CUdevice device;
    cuDeviceGet(&device, 0);
    CUcontext context;
    cuDevicePrimaryCtxRetain(&context, device);
    cuCtxPushCurrent(context);

    OptixDeviceContext optixDeviceContext = InitializeOptix(context);
    BVH *bvh;
    OptixAccelBuildOptions options = {};

    // options.operation
    //
    //     OPTIX_ASSERT(optixAccelBuild(
    //         optixDeviceContext, 0, &options, const OptixBuildInput *buildInputs,
    //         unsigned int numBuildInputs, CUdeviceptr tempBuffer, size_t
    //         tempBufferSizeInBytes, CUdeviceptr outputBuffer, size_t outputBufferSizeInBytes,
    //         OptixTraversableHandle *outputHandle, const OptixAccelEmitDesc
    //         *emittedProperties, unsigned int numEmittedProperties));

#if (OPTIX_VERSION >= 90000)
    // if (bvh->flags & BVHFlags::USE_CLUSTERS)
    {
        ClusterAccelerationStructureLimits limits;
        optixDeviceContextGetProperty(optixDeviceContext,
                                      OPTIX_DEVICE_PROPERTY_LIMIT_MAX_CLUSTER_TRIANGLES,
                                      &limits.maxTrianglesPerCluster, sizeof(unsigned int));
        optixDeviceContextGetProperty(optixDeviceContext,
                                      OPTIX_DEVICE_PROPERTY_LIMIT_MAX_CLUSTER_VERTICES,
                                      &limits.maxVerticesPerCluster, sizeof(unsigned int));
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

        optixClusterAccelBuild(optixDeviceContext, 0, &buildModeDesc, &buildInput,
                               CUdeviceptr(args), CUdeviceptr(argsCount), 0);
    }
#endif

    cuCtxPopCurrent(&context);
}

#endif

YBI_NAMESPACE_END
