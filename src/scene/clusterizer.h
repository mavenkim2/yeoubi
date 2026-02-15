#pragma once

#include "scene/scene.h"
#include <meshoptimizer.h>
#include <vector>

YBI_NAMESPACE_BEGIN

struct MeshletClustersResult
{
    size_t meshletCount;
    std::vector<meshopt_Meshlet> meshlets;
    std::vector<unsigned int> meshletVertices;
    std::vector<unsigned char> meshletTriangles;
};

inline MeshletClustersResult Test(const Mesh &mesh)
{
    const size_t indexCount = mesh.indices.size();
    const size_t vertexCount = mesh.positions.size();
    if (indexCount == 0 || vertexCount == 0)
    {
        return {};
    }

    const size_t maxVertices = 256;
    const size_t minTriangles = 1;
    const size_t maxTriangles = 128;
    const float fillWeight = 0.5f;

    const size_t meshletCountBound =
        meshopt_buildMeshletsBound(indexCount, maxVertices, minTriangles);

    MeshletClustersResult result;
    result.meshlets.resize(meshletCountBound);
    result.meshletVertices.resize(indexCount);
    result.meshletTriangles.resize(indexCount);

    const float *vertexPositions = reinterpret_cast<const float *>(mesh.positions.data());
    const size_t vertexPositionsStride = sizeof(float3);

    result.meshletCount = meshopt_buildMeshletsSpatial(result.meshlets.data(),
                                                       result.meshletVertices.data(),
                                                       result.meshletTriangles.data(),
                                                       mesh.indices.data(),
                                                       indexCount,
                                                       vertexPositions,
                                                       vertexCount,
                                                       vertexPositionsStride,
                                                       maxVertices,
                                                       minTriangles,
                                                       maxTriangles,
                                                       fillWeight);

    result.meshlets.resize(result.meshletCount);
    return result;
}

YBI_NAMESPACE_END
