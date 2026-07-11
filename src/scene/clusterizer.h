#pragma once

#include "core/assert.h"
#include "scene/scene.h"
#include <vector>

namespace ybi
{

struct DataContainer;
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
                    MeshletClustersResult &result);
void BuildMeshletIndexBuffer(const MeshletClustersResult &clusters, Array<int> &indices);

} // namespace ybi
