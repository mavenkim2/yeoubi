#pragma once

#include "core/assert.h"
#include "scene/scene.h"
#include <vector>

namespace ybi
{

struct DataContainer;
struct DeviceCluster;
class Mesh;

struct Meshlet
{
    uint32_t vertexOffset;
    uint32_t triangleOffset;

    uint32_t vertexCount;
    uint32_t triangleCount;
};

struct MeshletClustersResult
{
    Array<Meshlet> meshlets;
    Array<unsigned int> meshletVertices;
    Array<unsigned char> meshletTriangles;
};

void ClusterizeMesh(const DataContainer *container,
                    const Mesh &mesh,
                    const Array<int> &meshIndices,
                    MeshletClustersResult &result);
void BuildMeshletTriangleBuffer(const MeshletClustersResult &clusters,
                                const Array<int> &sourceIndices,
                                Array<int> &triangles);
void CompressClusters(MeshletClustersResult &clusters,
                      const Mesh &mesh,
                      DataContainer *container,
                      Array<DeviceCluster> &deviceClusters,
                      size_t clusterOffset);
void Help();

} // namespace ybi
