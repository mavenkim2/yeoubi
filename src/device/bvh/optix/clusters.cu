#include "device/bvh/optix/clusters.cuh"
#include "device/cuda_assert.h"
#include <cuda.h>

YBI_NAMESPACE_BEGIN

__global__ void
WriteTriangleClusterDescriptors(const MeshletDesc *meshlets,
                                const unsigned int *meshletVertices,
                                const uint8_t *meshletTriangles,
                                const float3 *sourcePositions,
                                unsigned int numClusters,
                                unsigned int *nextVertexOffset,
                                unsigned int *nextIndexOffset,
                                float3 *outputVertices,
                                uint8_t *outputIndices,
                                OptixClusterAccelBuildInputTrianglesArgs *buildArgs)
{
    for (uint32_t clusterIdx = blockIdx.x; clusterIdx < numClusters; clusterIdx += gridDim.x)
    {
        const MeshletDesc meshlet = meshlets[clusterIdx];
        const unsigned int vertexCount = meshlet.vertexCount;
        const unsigned int triangleCount = meshlet.triangleCount;

        const unsigned int vertexOffset = atomicAdd(nextVertexOffset, vertexCount);
        const unsigned int indexOffset = atomicAdd(nextIndexOffset, triangleCount * 3u);

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
    }
}

YBI_NAMESPACE_END
