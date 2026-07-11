#include "scene/clusterizer.h"
#include "core/thread_context.h"
#include <meshoptimizer.h>

namespace ybi
{

void ClusterizeMesh(const DataContainer *container,
                    const Mesh &mesh,
                    MeshletClustersResult &result)
{
    const Array<Vec3f> &meshPositions = container->GetMeshPositions(mesh);
    const Array<int> &meshIndices = container->GetMeshIndices(mesh);
    const size_t indexCount = meshIndices.size();
    const size_t vertexCount = meshPositions.size();

    YBI_ASSERT(indexCount);
    YBI_ASSERT(vertexCount);

    const size_t maxVertices = 256;
    const size_t minTriangles = 64;
    const size_t maxTriangles = 128;
    const float fillWeight = 0.5f;

    const size_t meshletCountBound =
        meshopt_buildMeshletsBound(indexCount, maxVertices, minTriangles);

    ScratchArena scratch;

    View<meshopt_Meshlet> tempMeshlets =
        scratch.arena->PushArray<meshopt_Meshlet>(meshletCountBound);

    result.meshletVertices.Resize(indexCount);
    result.meshletTriangles.Resize(indexCount);

    const float *vertexPositions = reinterpret_cast<const float *>(meshPositions.data());
    const size_t vertexPositionsStride = sizeof(meshPositions[0]);

    const size_t meshletCount = meshopt_buildMeshletsSpatial(tempMeshlets.data(),
                                                             result.meshletVertices.data(),
                                                             result.meshletTriangles.data(),
                                                             meshIndices.data(),
                                                             indexCount,
                                                             vertexPositions,
                                                             vertexCount,
                                                             vertexPositionsStride,
                                                             maxVertices,
                                                             minTriangles,
                                                             maxTriangles,
                                                             fillWeight);

    result.meshlets.Resize(meshletCount);
    for (size_t i = 0; i < meshletCount; i++)
    {
        Meshlet &meshlet = result.meshlets[i];
        meshopt_Meshlet &tempMeshlet = tempMeshlets[i];

        meshlet.vertexOffset = tempMeshlet.vertex_offset;
        meshlet.triangleOffset = tempMeshlet.triangle_offset;
        meshlet.vertexCount = tempMeshlet.vertex_count;
        meshlet.triangleCount = tempMeshlet.triangle_count;
    }
}

void BuildMeshletIndexBuffer(const MeshletClustersResult &clusters, Array<int> &indices)
{
    size_t indexCount = 0;
    for (size_t meshletIndex = 0; meshletIndex < clusters.meshlets.size(); meshletIndex++)
    {
        const Meshlet &meshlet = clusters.meshlets[meshletIndex];
        YBI_ASSERT(meshlet.vertexOffset + meshlet.vertexCount <= clusters.meshletVertices.size());
        YBI_ASSERT(meshlet.triangleOffset + meshlet.triangleCount * 3 <=
                   clusters.meshletTriangles.size());
        indexCount += meshlet.triangleCount * 3;
    }

    indices.Resize(indexCount);

    size_t outIndex = 0;
    for (size_t meshletIndex = 0; meshletIndex < clusters.meshlets.size(); meshletIndex++)
    {
        const Meshlet &meshlet = clusters.meshlets[meshletIndex];
        const unsigned int *meshletVertices =
            clusters.meshletVertices.data() + meshlet.vertexOffset;
        const unsigned char *meshletTriangles =
            clusters.meshletTriangles.data() + meshlet.triangleOffset;

        for (uint32_t triangleIndex = 0; triangleIndex < meshlet.triangleCount; triangleIndex++)
        {
            for (uint32_t cornerIndex = 0; cornerIndex < 3; cornerIndex++)
            {
                const unsigned char localVertexIndex = meshletTriangles[triangleIndex * 3 +
                                                                        cornerIndex];
                YBI_ASSERT(localVertexIndex < meshlet.vertexCount);
                indices[outIndex++] = static_cast<int>(meshletVertices[localVertexIndex]);
            }
        }
    }

    YBI_ASSERT(outIndex == indices.size());
}

} // namespace ybi
