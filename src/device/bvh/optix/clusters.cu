#include "device/bvh/optix/clusters.cuh"
#include "device/cuda_assert.h"
#include "device/cuda_device.h"
#include "optix_types.h"
#include "scene/clusterizer.h"
#include <algorithm>
#include <cuda.h>

YBI_NAMESPACE_BEGIN

void BuildMeshCLAS(CUDADevice *cudaDevice,
                   HostMemoryArena &hostArena,
                   const Mesh &mesh,
                   const MeshletClustersResult &clusterResult,
                   size_t &totalOutputSizeOut)
{
    const size_t meshletCount = clusterResult.meshletCount;
    const uint32_t vertexCount = static_cast<uint32_t>(mesh.positions.size());

    size_t totalVertexBound = 0;
    size_t totalIndexBound = 0;
    for (size_t i = 0; i < meshletCount; i++)
    {
        totalVertexBound += clusterResult.meshlets[i].vertex_count;
        totalIndexBound += clusterResult.meshlets[i].triangle_count * 3u;
    }

    CUDAMemoryArena &deviceArena = *cudaDevice->deviceArena;

    DeviceMemoryView<::float3> devicePositions = deviceArena.PushArray<::float3>(vertexCount);
    CUDA_ASSERT(cuMemcpyHtoD(CUdeviceptr(devicePositions.data()),
                             mesh.positions.data(),
                             vertexCount * sizeof(::float3)));

    MemoryView<MeshletDesc> hostMeshletDescs = hostArena.PushArray<MeshletDesc>(meshletCount);
    for (size_t i = 0; i < meshletCount; i++)
    {
        hostMeshletDescs[i].vertexOffset = clusterResult.meshlets[i].vertex_offset;
        hostMeshletDescs[i].triangleOffset = clusterResult.meshlets[i].triangle_offset;
        hostMeshletDescs[i].vertexCount = clusterResult.meshlets[i].vertex_count;
        hostMeshletDescs[i].triangleCount = clusterResult.meshlets[i].triangle_count;
    }

    const size_t meshletVerticesSize = clusterResult.meshletVertices.size();
    const size_t meshletTrianglesSize = clusterResult.meshletTriangles.size();

    DeviceMemoryView<MeshletDesc> deviceMeshlets =
        deviceArena.PushArray<MeshletDesc>(meshletCount);
    DeviceMemoryView<unsigned int> deviceMeshletVertices =
        deviceArena.PushArray<unsigned int>(meshletVerticesSize);
    DeviceMemoryView<uint8_t> deviceMeshletTriangles =
        deviceArena.PushArray<uint8_t>(meshletTrianglesSize);
    DeviceMemoryView<unsigned int> nextVertexOffset = deviceArena.PushArray<unsigned int>(1);
    DeviceMemoryView<unsigned int> nextIndexOffset = deviceArena.PushArray<unsigned int>(1);

    CUDA_ASSERT(cuMemcpyHtoD(CUdeviceptr(deviceMeshlets.data()),
                             hostMeshletDescs.data(),
                             meshletCount * sizeof(MeshletDesc)));
    CUDA_ASSERT(cuMemcpyHtoD(CUdeviceptr(deviceMeshletVertices.data()),
                             clusterResult.meshletVertices.data(),
                             meshletVerticesSize * sizeof(unsigned int)));
    CUDA_ASSERT(cuMemcpyHtoD(CUdeviceptr(deviceMeshletTriangles.data()),
                             clusterResult.meshletTriangles.data(),
                             meshletTrianglesSize));

    CUDA_ASSERT(cuMemsetD32(CUdeviceptr(nextVertexOffset.data()), 0, 1));
    CUDA_ASSERT(cuMemsetD32(CUdeviceptr(nextIndexOffset.data()), 0, 1));

    DeviceMemoryView<::float3> outputVertices = deviceArena.PushArray<::float3>(totalVertexBound);
    DeviceMemoryView<uint8_t> outputIndices = deviceArena.PushArray<uint8_t>(totalIndexBound);
    DeviceMemoryView<OptixClusterAccelBuildInputTrianglesArgs> buildArgs =
        deviceArena.PushArray<OptixClusterAccelBuildInputTrianglesArgs>(meshletCount);

    WriteTriangleClusterDescriptors<<<static_cast<unsigned int>(meshletCount), 256>>>(
        deviceMeshlets.data(),
        deviceMeshletVertices.data(),
        deviceMeshletTriangles.data(),
        devicePositions.data(),
        static_cast<unsigned int>(meshletCount),
        nextVertexOffset.data(),
        nextIndexOffset.data(),
        outputVertices.data(),
        outputIndices.data(),
        buildArgs.data());

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
    DeviceMemoryView<uint8_t> temp =
        deviceArena.PushArray<uint8_t>(tempSizeInBytes, OPTIX_ACCEL_BUFFER_BYTE_ALIGNMENT);
    DeviceMemoryView<uint32_t> deviceOutputSizes = deviceArena.PushArray<uint32_t>(meshletCount);

    DeviceMemoryView<uint32_t> deviceArgsCount = deviceArena.PushArray<uint32_t>(1);
    const uint32_t argsCountVal = static_cast<unsigned int>(meshletCount);
    CUDA_ASSERT(
        cuMemcpyHtoD(CUdeviceptr(deviceArgsCount.data()), &argsCountVal, sizeof(uint32_t)));

    OptixClusterAccelBuildModeDesc getSizeDesc = {};
    getSizeDesc.mode = OPTIX_CLUSTER_ACCEL_BUILD_MODE_GET_SIZES;
    getSizeDesc.getSize.outputSizesBuffer = CUdeviceptr(deviceOutputSizes.data());
    getSizeDesc.getSize.outputSizesStrideInBytes = sizeof(uint32_t);
    getSizeDesc.getSize.tempBuffer = CUdeviceptr(temp.data());
    getSizeDesc.getSize.tempBufferSizeInBytes = getSizeBufferSizes.tempSizeInBytes;

    OPTIX_ASSERT(optixClusterAccelBuild(cudaDevice->optixDeviceContext,
                                        0,
                                        &getSizeDesc,
                                        &buildInput,
                                        CUdeviceptr(buildArgs.data()),
                                        CUdeviceptr(deviceArgsCount.data()),
                                        sizeof(OptixClusterAccelBuildInputTrianglesArgs)));

    DeviceMemoryView<uint32_t> totalCLASSize = deviceArena.PushArray<uint32_t>(1);
    ComputeCLASTotalSize<<<1, 1>>>(
        deviceOutputSizes.data(), static_cast<uint32_t>(meshletCount), totalCLASSize.data());

    uint32_t totalSize = 0;
    CUDA_ASSERT(cuMemcpyDtoH(&totalSize, CUdeviceptr(totalCLASSize.data()), sizeof(uint32_t)));
    CUDA_ASSERT(cuStreamSynchronize(0));

    YBI_ASSERT(totalSize);

    DeviceMemoryView<uint8_t> outputBuffer = cudaDevice->Alloc<uint8_t>(totalSize);
    const CUdeviceptr baseAddr = CUdeviceptr(outputBuffer.data());

    DeviceMemoryView<uint64_t> deviceDestAddresses = deviceArena.PushArray<uint64_t>(meshletCount);
    ComputeCLASDestAddresses<<<1, 1>>>(deviceOutputSizes.data(),
                                       static_cast<uint32_t>(meshletCount),
                                       baseAddr,
                                       deviceDestAddresses.data());

    OptixClusterAccelBuildModeDesc explicitDesc = {};
    explicitDesc.mode = OPTIX_CLUSTER_ACCEL_BUILD_MODE_EXPLICIT_DESTINATIONS;
    explicitDesc.explicitDest.tempBuffer = CUdeviceptr(temp.data());
    explicitDesc.explicitDest.tempBufferSizeInBytes = explicitBufferSizes.tempSizeInBytes;
    explicitDesc.explicitDest.outputSizesBuffer = CUdeviceptr(deviceOutputSizes.data());
    explicitDesc.explicitDest.destAddressesBuffer = CUdeviceptr(deviceDestAddresses.data());
    explicitDesc.explicitDest.destAddressesStrideInBytes = sizeof(uint64_t);
    explicitDesc.explicitDest.outputHandlesBuffer = CUdeviceptr(deviceDestAddresses.data());

    OPTIX_ASSERT(optixClusterAccelBuild(cudaDevice->optixDeviceContext,
                                        0,
                                        &explicitDesc,
                                        &buildInput,
                                        CUdeviceptr(buildArgs.data()),
                                        CUdeviceptr(deviceArgsCount.data()),
                                        sizeof(OptixClusterAccelBuildInputTrianglesArgs)));
    CUDA_ASSERT(cuStreamSynchronize(0));

    // GAS BUILD
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
    DeviceMemoryView<uint8_t> gasTemp =
        deviceArena.PushArray<uint8_t>(gasTempSizeInBytes, OPTIX_ACCEL_BUFFER_BYTE_ALIGNMENT);
    DeviceMemoryView<uint32_t> deviceGASOutputSizes = deviceArena.PushArray<uint32_t>(1);

    DeviceMemoryView<OptixClusterAccelBuildInputClustersArgs> deviceClusterArgs =
        deviceArena.PushArray<OptixClusterAccelBuildInputClustersArgs>(1);
    OptixClusterAccelBuildInputClustersArgs hostClusterArgs = {};
    hostClusterArgs.clusterHandlesBuffer = CUdeviceptr(deviceDestAddresses.data());
    hostClusterArgs.clusterHandlesCount = static_cast<unsigned int>(meshletCount);
    hostClusterArgs.clusterHandlesBufferStrideInBytes = sizeof(OptixTraversableHandle);
    CUDA_ASSERT(cuMemcpyHtoD(CUdeviceptr(deviceClusterArgs.data()),
                             &hostClusterArgs,
                             sizeof(OptixClusterAccelBuildInputClustersArgs)));

    DeviceMemoryView<uint32_t> deviceGASArgsCount = deviceArena.PushArray<uint32_t>(1);
    const uint32_t gasArgsCountVal = 1u;
    CUDA_ASSERT(
        cuMemcpyHtoD(CUdeviceptr(deviceGASArgsCount.data()), &gasArgsCountVal, sizeof(uint32_t)));

    OptixClusterAccelBuildModeDesc gasGetSizeDesc = {};
    gasGetSizeDesc.mode = OPTIX_CLUSTER_ACCEL_BUILD_MODE_GET_SIZES;
    gasGetSizeDesc.getSize.outputSizesBuffer = CUdeviceptr(deviceGASOutputSizes.data());
    gasGetSizeDesc.getSize.outputSizesStrideInBytes = sizeof(uint32_t);
    gasGetSizeDesc.getSize.tempBuffer = CUdeviceptr(gasTemp.data());
    gasGetSizeDesc.getSize.tempBufferSizeInBytes = gasGetSizeBufferSizes.tempSizeInBytes;

    OPTIX_ASSERT(optixClusterAccelBuild(cudaDevice->optixDeviceContext,
                                        0,
                                        &gasGetSizeDesc,
                                        &gasBuildInput,
                                        CUdeviceptr(deviceClusterArgs.data()),
                                        CUdeviceptr(deviceGASArgsCount.data()),
                                        sizeof(OptixClusterAccelBuildInputClustersArgs)));

    uint32_t gasSize = 0;
    CUDA_ASSERT(
        cuMemcpyDtoH(&gasSize, CUdeviceptr(deviceGASOutputSizes.data()), sizeof(uint32_t)));

    YBI_ASSERT(gasSize);

    DeviceMemoryView<uint8_t> gasOutputBuffer = cudaDevice->Alloc<uint8_t>(gasSize);
    const uint64_t gasBaseAddr = reinterpret_cast<uint64_t>(gasOutputBuffer.data());
    DeviceMemoryView<uint64_t> deviceGASDestAddresses = deviceArena.PushArray<uint64_t>(1);
    CUDA_ASSERT(
        cuMemcpyHtoD(CUdeviceptr(deviceGASDestAddresses.data()), &gasBaseAddr, sizeof(uint64_t)));

    OptixClusterAccelBuildModeDesc gasExplicitDesc = {};
    gasExplicitDesc.mode = OPTIX_CLUSTER_ACCEL_BUILD_MODE_EXPLICIT_DESTINATIONS;
    gasExplicitDesc.explicitDest.tempBuffer = CUdeviceptr(gasTemp.data());
    gasExplicitDesc.explicitDest.tempBufferSizeInBytes = gasExplicitBufferSizes.tempSizeInBytes;
    gasExplicitDesc.explicitDest.outputSizesBuffer = CUdeviceptr(deviceGASOutputSizes.data());
    gasExplicitDesc.explicitDest.destAddressesBuffer = CUdeviceptr(deviceGASDestAddresses.data());
    gasExplicitDesc.explicitDest.destAddressesStrideInBytes = sizeof(uint64_t);

    OPTIX_ASSERT(optixClusterAccelBuild(cudaDevice->optixDeviceContext,
                                        0,
                                        &gasExplicitDesc,
                                        &gasBuildInput,
                                        CUdeviceptr(deviceClusterArgs.data()),
                                        CUdeviceptr(deviceGASArgsCount.data()),
                                        sizeof(OptixClusterAccelBuildInputClustersArgs)));
    CUDA_ASSERT(cuStreamSynchronize(0));

    cudaDevice->bvhTotalAllocated += gasSize;
    cudaDevice->bvhTotalAllocated += totalSize;
    totalOutputSizeOut = totalSize;

    deviceArena.Clear();
}

__global__ void
WriteTriangleClusterDescriptors(const MeshletDesc *meshlets,
                                const unsigned int *meshletVertices,
                                const uint8_t *meshletTriangles,
                                const ::float3 *sourcePositions,
                                unsigned int numClusters,
                                unsigned int *nextVertexOffset,
                                unsigned int *nextIndexOffset,
                                ::float3 *outputVertices,
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
            args.vertexBufferStrideInBytes = 0;
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
        return;

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
        return;

    uint64_t offset = 0;
    for (uint32_t i = 0; i < numClusters; i++)
    {
        destAddresses[i] = baseAddr + offset;
        offset += sizes[i];
    }
}

YBI_NAMESPACE_END
