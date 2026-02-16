#pragma once

#include "util/host_memory_arena.h"
#include <cstdint>
#include <cuda.h>
#include <cuda_runtime.h>
#include <optix_types.h>

YBI_NAMESPACE_BEGIN

struct CUDADevice;
struct MeshletClustersResult;
struct Mesh;

void BuildMeshCLAS(CUDADevice *cudaDevice,
                   HostMemoryArena &hostArena,
                   const Mesh &mesh,
                   const MeshletClustersResult &clusterResult,
                   size_t &totalOutputSizeOut);

struct MeshletDesc
{
    unsigned int vertexOffset;
    unsigned int triangleOffset;
    unsigned int vertexCount;
    unsigned int triangleCount;
};

extern "C" __global__ void
ComputeCLASTotalSize(const uint32_t *sizes, uint32_t numClusters, uint32_t *totalSizeOut);

extern "C" __global__ void ComputeCLASDestAddresses(const uint32_t *sizes,
                                                    uint32_t numClusters,
                                                    uint64_t baseAddr,
                                                    uint64_t *destAddresses);

extern "C" __global__ void
WriteTriangleClusterDescriptors(const MeshletDesc *meshlets,
                                const unsigned int *meshletVertices,
                                const uint8_t *meshletTriangles,
                                const float3 *sourcePositions,
                                unsigned int numClusters,
                                unsigned int *nextVertexOffset,
                                unsigned int *nextIndexOffset,
                                float3 *outputVertices,
                                uint8_t *outputIndices,
                                OptixClusterAccelBuildInputTrianglesArgs *buildArgs);

YBI_NAMESPACE_END
