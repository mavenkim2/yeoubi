#include "core/assert.h"
#include "core/vec3f.h"
#include "device/cuda_assert.h"
#include "device/cuda_device.h"
#include "scene/clusterizer.h"
#include "scene/scene.h"
#include <algorithm>
#include <cstdint>
#include <cuda.h>
#include <optix_types.h>
#include <vector>

namespace ybi
{

__global__ void
WriteTriangleClusterDescriptors(const Meshlet *meshlets,
                                const unsigned int *meshletVertices,
                                const uint8_t *meshletTriangles,
                                const Vec3f *sourcePositions,
                                unsigned int numClusters,
                                unsigned int *nextVertexOffset,
                                unsigned int *nextIndexOffset,
                                Vec3f *outputVertices,
                                uint8_t *outputIndices,
                                OptixClusterAccelBuildInputTrianglesArgs *buildArgs);

__global__ void
ComputeCLASTotalSize(const uint32_t *sizes, uint32_t numClusters, uint32_t *totalSizeOut);

__global__ void ComputeCLASDestAddresses(const uint32_t *sizes,
                                         uint32_t numClusters,
                                         uint64_t baseAddr,
                                         uint64_t *destAddresses);

void BuildMeshCLAS(CUDADevice *cudaDevice,
                   DataContainer *container,
                   const Mesh &mesh,
                   const MeshletClustersResult &clusterResult,
                   size_t &totalOutputSizeOut,
                   OptixTraversableHandle &outputHandleOut)
{
    YBI_ASSERT(cudaDevice);
    YBI_ASSERT(container);

    totalOutputSizeOut = 0;
    outputHandleOut = {};

    const BufferHandle positionHandle = mesh.GetPositionsHandle();
    YBI_ASSERT(positionHandle.IsValid());
    YBI_ASSERT(positionHandle.index < container->meshPositions.size());

    const Array<Vec3f> &positions = container->meshPositions[positionHandle.index];
    const size_t meshletCount = clusterResult.meshletCount;
    const uint32_t vertexCount = static_cast<uint32_t>(positions.size());
    if (meshletCount == 0 || vertexCount == 0)
    {
        return;
    }

    size_t totalVertexBound = 0;
    size_t totalIndexBound = 0;
    for (size_t i = 0; i < meshletCount; i++)
    {
        totalVertexBound += clusterResult.meshlets[i].vertex_count;
        totalIndexBound += clusterResult.meshlets[i].triangle_count * 3u;
    }

    std::vector<DeviceBuffer> temporaryBuffers;
    auto allocTemporary = [&](size_t bytes) {
        DeviceBuffer buffer = cudaDevice->Alloc(bytes);
        temporaryBuffers.push_back(buffer);
        return buffer;
    };

    DeviceBuffer devicePositions = allocTemporary(vertexCount * sizeof(Vec3f));
    CUDA_ASSERT(cuMemcpyHtoD(
        CUdeviceptr(devicePositions.Data()), positions.data(), vertexCount * sizeof(Vec3f)));

    std::vector<Meshlet> hostMeshletDescs(meshletCount);
    for (size_t i = 0; i < meshletCount; i++)
    {
        hostMeshletDescs[i].vertexOffset = clusterResult.meshlets[i].vertex_offset;
        hostMeshletDescs[i].triangleOffset = clusterResult.meshlets[i].triangle_offset;
        hostMeshletDescs[i].vertexCount = clusterResult.meshlets[i].vertex_count;
        hostMeshletDescs[i].triangleCount = clusterResult.meshlets[i].triangle_count;
    }

    const size_t meshletVerticesSize = clusterResult.meshletVertices.size();
    const size_t meshletTrianglesSize = clusterResult.meshletTriangles.size();

    DeviceBuffer deviceMeshlets = allocTemporary(meshletCount * sizeof(Meshlet));
    DeviceBuffer deviceMeshletVertices =
        allocTemporary(meshletVerticesSize * sizeof(unsigned int));
    DeviceBuffer deviceMeshletTriangles = allocTemporary(meshletTrianglesSize * sizeof(uint8_t));
    DeviceBuffer nextVertexOffset = allocTemporary(sizeof(unsigned int));
    DeviceBuffer nextIndexOffset = allocTemporary(sizeof(unsigned int));

    CUDA_ASSERT(cuMemcpyHtoD(CUdeviceptr(deviceMeshlets.Data()),
                             hostMeshletDescs.data(),
                             meshletCount * sizeof(Meshlet)));
    CUDA_ASSERT(cuMemcpyHtoD(CUdeviceptr(deviceMeshletVertices.Data()),
                             clusterResult.meshletVertices.data(),
                             meshletVerticesSize * sizeof(unsigned int)));
    CUDA_ASSERT(cuMemcpyHtoD(CUdeviceptr(deviceMeshletTriangles.Data()),
                             clusterResult.meshletTriangles.data(),
                             meshletTrianglesSize));

    CUDA_ASSERT(cuMemsetD32(CUdeviceptr(nextVertexOffset.Data()), 0, 1));
    CUDA_ASSERT(cuMemsetD32(CUdeviceptr(nextIndexOffset.Data()), 0, 1));

    DeviceBuffer outputVertices = allocTemporary(totalVertexBound * sizeof(Vec3f));
    DeviceBuffer outputIndices = allocTemporary(totalIndexBound * sizeof(uint8_t));
    DeviceBuffer buildArgs =
        allocTemporary(meshletCount * sizeof(OptixClusterAccelBuildInputTrianglesArgs));

    WriteTriangleClusterDescriptors<<<static_cast<unsigned int>(meshletCount), 256>>>(
        deviceMeshlets.As<MeshletDesc>().Get(),
        deviceMeshletVertices.As<unsigned int>().Get(),
        deviceMeshletTriangles.As<uint8_t>().Get(),
        devicePositions.As<Vec3f>().Get(),
        static_cast<unsigned int>(meshletCount),
        nextVertexOffset.As<unsigned int>().Get(),
        nextIndexOffset.As<unsigned int>().Get(),
        outputVertices.As<Vec3f>().Get(),
        outputIndices.As<uint8_t>().Get(),
        buildArgs.As<OptixClusterAccelBuildInputTrianglesArgs>().Get());

    OptixClusterAccelBuildInput buildInput = {};
    buildInput.type = OPTIX_CLUSTER_ACCEL_BUILD_TYPE_CLUSTERS_FROM_TRIANGLES;
    buildInput.triangles.flags = OPTIX_CLUSTER_ACCEL_BUILD_FLAG_PREFER_FAST_TRACE;
    buildInput.triangles.maxArgCount = static_cast<unsigned int>(meshletCount);
    buildInput.triangles.vertexFormat = OPTIX_VERTEX_FORMAT_FLOAT3;
    buildInput.triangles.maxSbtIndexValue = 0;
    buildInput.triangles.maxUniqueSbtIndexCountPerArg = 1;
    buildInput.triangles.maxTriangleCountPerArg = 128;
    buildInput.triangles.maxVertexCountPerArg = 256;

    OptixAccelBufferSizes getSizeBufferSizes = {};
    OPTIX_ASSERT(optixClusterAccelComputeMemoryUsage(cudaDevice->optixDeviceContext,
                                                     OPTIX_CLUSTER_ACCEL_BUILD_MODE_GET_SIZES,
                                                     &buildInput,
                                                     &getSizeBufferSizes));

    OptixAccelBufferSizes explicitBufferSizes = {};
    OPTIX_ASSERT(
        optixClusterAccelComputeMemoryUsage(cudaDevice->optixDeviceContext,
                                            OPTIX_CLUSTER_ACCEL_BUILD_MODE_EXPLICIT_DESTINATIONS,
                                            &buildInput,
                                            &explicitBufferSizes));

    const size_t tempSizeInBytes =
        std::max(getSizeBufferSizes.tempSizeInBytes, explicitBufferSizes.tempSizeInBytes);
    DeviceBuffer temp = allocTemporary(tempSizeInBytes);
    DeviceBuffer deviceOutputSizes = allocTemporary(meshletCount * sizeof(uint32_t));

    DeviceBuffer deviceArgsCount = allocTemporary(sizeof(uint32_t));
    const uint32_t argsCountVal = static_cast<unsigned int>(meshletCount);
    CUDA_ASSERT(
        cuMemcpyHtoD(CUdeviceptr(deviceArgsCount.Data()), &argsCountVal, sizeof(uint32_t)));

    OptixClusterAccelBuildModeDesc getSizeDesc = {};
    getSizeDesc.mode = OPTIX_CLUSTER_ACCEL_BUILD_MODE_GET_SIZES;
    getSizeDesc.getSize.outputSizesBuffer = CUdeviceptr(deviceOutputSizes.Data());
    getSizeDesc.getSize.outputSizesStrideInBytes = sizeof(uint32_t);
    getSizeDesc.getSize.tempBuffer = CUdeviceptr(temp.Data());
    getSizeDesc.getSize.tempBufferSizeInBytes = getSizeBufferSizes.tempSizeInBytes;

    OPTIX_ASSERT(optixClusterAccelBuild(cudaDevice->optixDeviceContext,
                                        0,
                                        &getSizeDesc,
                                        &buildInput,
                                        CUdeviceptr(buildArgs.Data()),
                                        CUdeviceptr(deviceArgsCount.Data()),
                                        sizeof(OptixClusterAccelBuildInputTrianglesArgs)));

    DeviceBuffer totalCLASSize = allocTemporary(sizeof(uint32_t));
    ComputeCLASTotalSize<<<1, 1>>>(deviceOutputSizes.As<uint32_t>().Get(),
                                   static_cast<uint32_t>(meshletCount),
                                   totalCLASSize.As<uint32_t>().Get());

    uint32_t totalSize = 0;
    CUDA_ASSERT(cuMemcpyDtoH(&totalSize, CUdeviceptr(totalCLASSize.Data()), sizeof(uint32_t)));
    CUDA_ASSERT(cuStreamSynchronize(0));

    YBI_ASSERT(totalSize);

    DeviceBuffer outputBuffer = cudaDevice->Alloc(totalSize);
    const CUdeviceptr baseAddr = CUdeviceptr(outputBuffer.Data());

    DeviceBuffer deviceDestAddresses = allocTemporary(meshletCount * sizeof(uint64_t));
    ComputeCLASDestAddresses<<<1, 1>>>(deviceOutputSizes.As<uint32_t>().Get(),
                                       static_cast<uint32_t>(meshletCount),
                                       baseAddr,
                                       deviceDestAddresses.As<uint64_t>().Get());

    OptixClusterAccelBuildModeDesc explicitDesc = {};
    explicitDesc.mode = OPTIX_CLUSTER_ACCEL_BUILD_MODE_EXPLICIT_DESTINATIONS;
    explicitDesc.explicitDest.tempBuffer = CUdeviceptr(temp.Data());
    explicitDesc.explicitDest.tempBufferSizeInBytes = explicitBufferSizes.tempSizeInBytes;
    explicitDesc.explicitDest.outputSizesBuffer = CUdeviceptr(deviceOutputSizes.Data());
    explicitDesc.explicitDest.destAddressesBuffer = CUdeviceptr(deviceDestAddresses.Data());
    explicitDesc.explicitDest.destAddressesStrideInBytes = sizeof(uint64_t);
    explicitDesc.explicitDest.outputHandlesBuffer = CUdeviceptr(deviceDestAddresses.Data());
    explicitDesc.explicitDest.outputHandlesStrideInBytes = sizeof(OptixTraversableHandle);

    OPTIX_ASSERT(optixClusterAccelBuild(cudaDevice->optixDeviceContext,
                                        0,
                                        &explicitDesc,
                                        &buildInput,
                                        CUdeviceptr(buildArgs.Data()),
                                        CUdeviceptr(deviceArgsCount.Data()),
                                        sizeof(OptixClusterAccelBuildInputTrianglesArgs)));
    CUDA_ASSERT(cuStreamSynchronize(0));

    OptixClusterAccelBuildInput gasBuildInput = {};
    gasBuildInput.type = OPTIX_CLUSTER_ACCEL_BUILD_TYPE_GASES_FROM_CLUSTERS;
    gasBuildInput.clusters.flags = OPTIX_CLUSTER_ACCEL_BUILD_FLAG_PREFER_FAST_TRACE;
    gasBuildInput.clusters.maxArgCount = 1;
    gasBuildInput.clusters.maxTotalClusterCount = static_cast<unsigned int>(meshletCount);
    gasBuildInput.clusters.maxClusterCountPerArg = static_cast<unsigned int>(meshletCount);

    OptixAccelBufferSizes gasGetSizeBufferSizes = {};
    OPTIX_ASSERT(optixClusterAccelComputeMemoryUsage(cudaDevice->optixDeviceContext,
                                                     OPTIX_CLUSTER_ACCEL_BUILD_MODE_GET_SIZES,
                                                     &gasBuildInput,
                                                     &gasGetSizeBufferSizes));

    OptixAccelBufferSizes gasExplicitBufferSizes = {};
    OPTIX_ASSERT(
        optixClusterAccelComputeMemoryUsage(cudaDevice->optixDeviceContext,
                                            OPTIX_CLUSTER_ACCEL_BUILD_MODE_EXPLICIT_DESTINATIONS,
                                            &gasBuildInput,
                                            &gasExplicitBufferSizes));

    const size_t gasTempSizeInBytes =
        std::max(gasGetSizeBufferSizes.tempSizeInBytes, gasExplicitBufferSizes.tempSizeInBytes);
    DeviceBuffer gasTemp = allocTemporary(gasTempSizeInBytes);
    DeviceBuffer deviceGASOutputSizes = allocTemporary(sizeof(uint32_t));

    DeviceBuffer deviceClusterArgs =
        allocTemporary(sizeof(OptixClusterAccelBuildInputClustersArgs));
    OptixClusterAccelBuildInputClustersArgs hostClusterArgs = {};
    hostClusterArgs.clusterHandlesBuffer = CUdeviceptr(deviceDestAddresses.Data());
    hostClusterArgs.clusterHandlesCount = static_cast<unsigned int>(meshletCount);
    hostClusterArgs.clusterHandlesBufferStrideInBytes = sizeof(OptixTraversableHandle);
    CUDA_ASSERT(cuMemcpyHtoD(CUdeviceptr(deviceClusterArgs.Data()),
                             &hostClusterArgs,
                             sizeof(OptixClusterAccelBuildInputClustersArgs)));

    DeviceBuffer deviceGASArgsCount = allocTemporary(sizeof(uint32_t));
    const uint32_t gasArgsCountVal = 1u;
    CUDA_ASSERT(
        cuMemcpyHtoD(CUdeviceptr(deviceGASArgsCount.Data()), &gasArgsCountVal, sizeof(uint32_t)));

    OptixClusterAccelBuildModeDesc gasGetSizeDesc = {};
    gasGetSizeDesc.mode = OPTIX_CLUSTER_ACCEL_BUILD_MODE_GET_SIZES;
    gasGetSizeDesc.getSize.outputSizesBuffer = CUdeviceptr(deviceGASOutputSizes.Data());
    gasGetSizeDesc.getSize.outputSizesStrideInBytes = sizeof(uint32_t);
    gasGetSizeDesc.getSize.tempBuffer = CUdeviceptr(gasTemp.Data());
    gasGetSizeDesc.getSize.tempBufferSizeInBytes = gasGetSizeBufferSizes.tempSizeInBytes;

    OPTIX_ASSERT(optixClusterAccelBuild(cudaDevice->optixDeviceContext,
                                        0,
                                        &gasGetSizeDesc,
                                        &gasBuildInput,
                                        CUdeviceptr(deviceClusterArgs.Data()),
                                        CUdeviceptr(deviceGASArgsCount.Data()),
                                        sizeof(OptixClusterAccelBuildInputClustersArgs)));

    uint32_t gasSize = 0;
    CUDA_ASSERT(
        cuMemcpyDtoH(&gasSize, CUdeviceptr(deviceGASOutputSizes.Data()), sizeof(uint32_t)));

    YBI_ASSERT(gasSize);

    DeviceBuffer gasOutputBuffer = cudaDevice->Alloc(gasSize);
    const uint64_t gasBaseAddr = reinterpret_cast<uint64_t>(gasOutputBuffer.Data());
    DeviceBuffer deviceGASDestAddresses = allocTemporary(sizeof(uint64_t));
    CUDA_ASSERT(
        cuMemcpyHtoD(CUdeviceptr(deviceGASDestAddresses.Data()), &gasBaseAddr, sizeof(uint64_t)));
    DeviceBuffer deviceGASOutputHandles = allocTemporary(sizeof(OptixTraversableHandle));

    OptixClusterAccelBuildModeDesc gasExplicitDesc = {};
    gasExplicitDesc.mode = OPTIX_CLUSTER_ACCEL_BUILD_MODE_EXPLICIT_DESTINATIONS;
    gasExplicitDesc.explicitDest.tempBuffer = CUdeviceptr(gasTemp.Data());
    gasExplicitDesc.explicitDest.tempBufferSizeInBytes = gasExplicitBufferSizes.tempSizeInBytes;
    gasExplicitDesc.explicitDest.outputSizesBuffer = CUdeviceptr(deviceGASOutputSizes.Data());
    gasExplicitDesc.explicitDest.destAddressesBuffer = CUdeviceptr(deviceGASDestAddresses.Data());
    gasExplicitDesc.explicitDest.destAddressesStrideInBytes = sizeof(uint64_t);
    gasExplicitDesc.explicitDest.outputHandlesBuffer = CUdeviceptr(deviceGASOutputHandles.Data());
    gasExplicitDesc.explicitDest.outputHandlesStrideInBytes = sizeof(OptixTraversableHandle);

    OPTIX_ASSERT(optixClusterAccelBuild(cudaDevice->optixDeviceContext,
                                        0,
                                        &gasExplicitDesc,
                                        &gasBuildInput,
                                        CUdeviceptr(deviceClusterArgs.Data()),
                                        CUdeviceptr(deviceGASArgsCount.Data()),
                                        sizeof(OptixClusterAccelBuildInputClustersArgs)));
    CUDA_ASSERT(cuStreamSynchronize(0));
    CUDA_ASSERT(cuMemcpyDtoH(&outputHandleOut,
                             CUdeviceptr(deviceGASOutputHandles.Data()),
                             sizeof(OptixTraversableHandle)));

    cudaDevice->bvhTotalAllocated += gasSize;
    cudaDevice->bvhTotalAllocated += totalSize;
    totalOutputSizeOut = totalSize;

    for (DeviceBuffer &buffer : temporaryBuffers)
    {
        cudaDevice->Free(buffer);
    }
}

__global__ void
WriteTriangleClusterDescriptors(const MeshletDesc *meshlets,
                                const unsigned int *meshletVertices,
                                const uint8_t *meshletTriangles,
                                const Vec3f *sourcePositions,
                                unsigned int numClusters,
                                unsigned int *nextVertexOffset,
                                unsigned int *nextIndexOffset,
                                Vec3f *outputVertices,
                                uint8_t *outputIndices,
                                OptixClusterAccelBuildInputTrianglesArgs *buildArgs)
{
    __shared__ unsigned int vertexOffset;
    __shared__ unsigned int indexOffset;

    for (uint32_t clusterIdx = blockIdx.x; clusterIdx < numClusters; clusterIdx += gridDim.x)
    {
        const MeshletDesc meshlet = meshlets[clusterIdx];
        const unsigned int vertexCount = meshlet.vertexCount;
        const unsigned int triangleCount = meshlet.triangleCount;

        if (threadIdx.x == 0)
        {
            vertexOffset = atomicAdd(nextVertexOffset, vertexCount);
            indexOffset = atomicAdd(nextIndexOffset, triangleCount * 3u);
        }
        __syncthreads();

        for (unsigned int i = threadIdx.x; i < vertexCount; i += blockDim.x)
        {
            const unsigned int globalVertexIdx = meshletVertices[meshlet.vertexOffset + i];
            outputVertices[vertexOffset + i] = sourcePositions[globalVertexIdx];
        }

        for (unsigned int t = threadIdx.x; t < triangleCount; t += blockDim.x)
        {
            const unsigned int src = meshlet.triangleOffset + t * 3u;
            outputIndices[indexOffset + t * 3u + 0u] = meshletTriangles[src + 0u];
            outputIndices[indexOffset + t * 3u + 1u] = meshletTriangles[src + 1u];
            outputIndices[indexOffset + t * 3u + 2u] = meshletTriangles[src + 2u];
        }

        if (threadIdx.x == 0)
        {
            OptixClusterAccelBuildInputTrianglesArgs &args = buildArgs[clusterIdx];
            args.clusterId = clusterIdx;
            args.clusterFlags = 0;
            args.triangleCount = triangleCount;
            args.vertexCount = vertexCount;
            args.positionTruncateBitCount = 0;
            args.indexFormat = OPTIX_CLUSTER_ACCEL_INDICES_FORMAT_8BIT;
            args.opacityMicromapIndexFormat = 0;
            args.basePrimitiveInfo = {};
            args.indexBufferStrideInBytes = 0;
            args.vertexBufferStrideInBytes = sizeof(Vec3f);
            args.primitiveInfoBufferStrideInBytes = 0;
            args.opacityMicromapIndexBufferStrideInBytes = 0;
            args.indexBuffer = reinterpret_cast<CUdeviceptr>(outputIndices + indexOffset);
            args.vertexBuffer = reinterpret_cast<CUdeviceptr>(outputVertices + vertexOffset);
            args.primitiveInfoBuffer = 0;
            args.opacityMicromapArray = 0;
            args.opacityMicromapIndexBuffer = 0;
            args.instantiationBoundingBoxLimit = 0;
        }
        __syncthreads();
    }
}

__global__ void
ComputeCLASTotalSize(const uint32_t *sizes, uint32_t numClusters, uint32_t *totalSizeOut)
{
    if (threadIdx.x > 0)
    {
        return;
    }

    uint32_t total = 0;
    for (uint32_t i = 0; i < numClusters; i++)
    {
        total += sizes[i];
    }
    *totalSizeOut = total;
}

__global__ void ComputeCLASDestAddresses(const uint32_t *sizes,
                                         uint32_t numClusters,
                                         uint64_t baseAddr,
                                         uint64_t *destAddresses)
{
    if (threadIdx.x > 0)
    {
        return;
    }

    uint64_t offset = 0;
    for (uint32_t i = 0; i < numClusters; i++)
    {
        destAddresses[i] = baseAddr + offset;
        offset += sizes[i];
    }
}

} // namespace ybi
