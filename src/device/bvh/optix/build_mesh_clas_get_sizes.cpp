#include "device/bvh/optix/build_mesh_clas_get_sizes.h"
#include "device/bvh/optix/clusters.cuh"
#include "device/cuda_device.h"

YBI_NAMESPACE_BEGIN

void BuildMeshCLASGetSizes(CUDADevice *cudaDevice,
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

    DeviceMemoryView<float3> devicePositions = deviceArena.PushArray<float3>(vertexCount);
    CUDA_ASSERT(cuMemcpyHtoD(
        CUdeviceptr(devicePositions.data()), mesh.positions.data(), vertexCount * sizeof(float3)));

    MemoryView<MeshletDesc> hostMeshletDescs = hostArena.PushArray<MeshletDesc>(meshletCount);
    for (size_t i = 0; i < meshletCount; i++)
    {
        hostMeshletDescs[i].vertexOffset = clusterResult.meshlets[i].vertex_offset;
        hostMeshletDescs[i].triangleOffset = clusterResult.meshlets[i].triangle_offset;
        hostMeshletDescs[i].vertexCount = clusterResult.meshlets[i].vertex_count;
        hostMeshletDescs[i].triangleCount = clusterResult.meshlets[i].triangle_count;
    }
    DeviceMemoryView<MeshletDesc> deviceMeshlets =
        deviceArena.PushArray<MeshletDesc>(meshletCount);
    CUDA_ASSERT(cuMemcpyHtoD(CUdeviceptr(deviceMeshlets.data()),
                             hostMeshletDescs.data(),
                             meshletCount * sizeof(MeshletDesc)));

    const size_t meshletVerticesSize = clusterResult.meshletVertices.size();
    DeviceMemoryView<unsigned int> deviceMeshletVertices =
        deviceArena.PushArray<unsigned int>(meshletVerticesSize);
    CUDA_ASSERT(cuMemcpyHtoD(CUdeviceptr(deviceMeshletVertices.data()),
                             clusterResult.meshletVertices.data(),
                             meshletVerticesSize * sizeof(unsigned int)));

    const size_t meshletTrianglesSize = clusterResult.meshletTriangles.size();
    DeviceMemoryView<uint8_t> deviceMeshletTriangles =
        deviceArena.PushArray<uint8_t>(meshletTrianglesSize);
    CUDA_ASSERT(cuMemcpyHtoD(CUdeviceptr(deviceMeshletTriangles.data()),
                             clusterResult.meshletTriangles.data(),
                             meshletTrianglesSize));

    DeviceMemoryView<unsigned int> nextVertexOffset = deviceArena.PushArray<unsigned int>(1);
    DeviceMemoryView<unsigned int> nextIndexOffset = deviceArena.PushArray<unsigned int>(1);
    CUDA_ASSERT(cuMemsetD32(CUdeviceptr(nextVertexOffset.data()), 0, 1));
    CUDA_ASSERT(cuMemsetD32(CUdeviceptr(nextIndexOffset.data()), 0, 1));

    DeviceMemoryView<float3> outputVertices = deviceArena.PushArray<float3>(totalVertexBound);
    DeviceMemoryView<uint8_t> outputIndices = deviceArena.PushArray<uint8_t>(totalIndexBound);
    DeviceMemoryView<OptixClusterAccelBuildInputTrianglesArgs> buildArgs =
        deviceArena.PushArray<OptixClusterAccelBuildInputTrianglesArgs>(meshletCount);

    WriteTriangleClusterDescriptors<<<meshletCount, 256>>>(deviceMeshlets.data(),
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

    DeviceMemoryView<uint8_t> temp = deviceArena.PushArray<uint8_t>(
        getSizeBufferSizes.tempSizeInBytes, OPTIX_ACCEL_BUFFER_BYTE_ALIGNMENT);
    DeviceMemoryView<uint32_t> deviceOutputSizes = deviceArena.PushArray<uint32_t>(meshletCount);

    DeviceMemoryView<uint32_t> deviceArgsCount = deviceArena.PushArray<uint32_t>(1);
    const uint32_t argsCountVal = static_cast<uint32_t>(meshletCount);
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
    CUDA_ASSERT(cuStreamSynchronize(0));

    MemoryView<uint32_t> hostOutputSizes = hostArena.PushArray<uint32_t>(meshletCount);
    CUDA_ASSERT(cuMemcpyDtoH(hostOutputSizes.data(),
                             CUdeviceptr(deviceOutputSizes.data()),
                             meshletCount * sizeof(uint32_t)));

    totalOutputSizeOut = 0;
    for (size_t i = 0; i < meshletCount; i++)
    {
        totalOutputSizeOut += hostOutputSizes[i];
    }
}

YBI_NAMESPACE_END
