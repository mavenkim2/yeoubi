#include "scene/clusterizer.h"
#include "core/assert.h"
#include "core/bit.h"
#include "core/bit_vector.h"
#include "core/math.h"
#include "core/octahedral.h"
#include "core/thread_context.h"
#include "core/vec3f.h"
#include "core/vec4u.h"
#include "kernel/geometry/cluster.h"
#include "scene/scene.h"
#include <atomic>
#include <cfloat>
#include <cstdint>
#include <functional>
#include <meshoptimizer.h>
#include <tbb/parallel_for.h>
#include <unordered_map>
#include <utility>

namespace ybi
{

#if defined(WITH_CUDA)
void ValidateCompressedClusterCUDA(const Array<uint8_t> &compressedBytes,
                                   const DeviceCluster &cluster,
                                   const DataContainer &container,
                                   const Mesh &mesh,
                                   const MeshletClustersResult &clusters,
                                   const Meshlet &meshlet,
                                   const Array<unsigned int> &vertexIndices,
                                   const Array<unsigned int> &faceVaryingIndices,
                                   const Array<int> &clusterTriangleToMeshTriangle,
                                   const Array<uint32_t> &triangleOrder,
                                   const Array<int> &triangleRotations,
                                   const Array<uint32_t> &oldIndexToNewIndex);
#endif

void ClusterizeMesh(const DataContainer *container,
                    const Mesh &mesh,
                    const Array<int> &meshIndices,
                    MeshletClustersResult &result)
{
    const Array<Vec3f> &meshPositions = container->GetMeshPositions(mesh);
    const size_t indexCount = meshIndices.size();
    const size_t vertexCount = meshPositions.size();

    YBI_ASSERT(indexCount);
    YBI_ASSERT(vertexCount);

    const size_t maxVertices = DEVICE_CLUSTER_MAX_VERTICES;
    const size_t minTriangles = DEVICE_CLUSTER_MAX_TRIANGLES / 2u;
    const size_t maxTriangles = DEVICE_CLUSTER_MAX_TRIANGLES;
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

void BuildMeshletTriangleBuffer(const MeshletClustersResult &clusters,
                                const Array<int> &sourceIndices,
                                Array<int> &triangles)
{
    struct TriangleKey
    {
        int indices[3];

        bool operator==(const TriangleKey &other) const
        {
            return indices[0] == other.indices[0] && indices[1] == other.indices[1] &&
                   indices[2] == other.indices[2];
        }
    };

    struct TriangleKeyHasher
    {
        size_t operator()(const TriangleKey &key) const
        {
            return static_cast<size_t>(MurmurHash64A(
                reinterpret_cast<const uint8_t *>(key.indices), sizeof(key.indices), 0));
        }
    };

    struct TriangleOccurrences
    {
        std::vector<int> sourceTriangles;
        size_t next = 0;
    };

    YBI_ASSERT(sourceIndices.size() % 3 == 0);

    const size_t sourceTriangleCount = sourceIndices.size() / 3;
    std::unordered_map<TriangleKey, TriangleOccurrences, TriangleKeyHasher> triangleMap;
    triangleMap.reserve(sourceTriangleCount);

    for (size_t triangleIndex = 0; triangleIndex < sourceTriangleCount; triangleIndex++)
    {
        const int *indices = sourceIndices.data() + triangleIndex * 3;
        const TriangleKey key = {{indices[0], indices[1], indices[2]}};
        triangleMap[key].sourceTriangles.push_back(static_cast<int>(triangleIndex));
    }

    size_t clusteredTriangleCount = 0;
    for (const Meshlet &meshlet : clusters.meshlets)
    {
        clusteredTriangleCount += meshlet.triangleCount;
    }
    YBI_ASSERT(clusteredTriangleCount <= sourceTriangleCount);
    triangles.Resize(clusteredTriangleCount);

    size_t outputTriangleIndex = 0;
    for (size_t meshletIndex = 0; meshletIndex < clusters.meshlets.size(); meshletIndex++)
    {
        const Meshlet &meshlet = clusters.meshlets[meshletIndex];
        YBI_ASSERT(meshlet.vertexOffset + meshlet.vertexCount <= clusters.meshletVertices.size());
        YBI_ASSERT(meshlet.triangleOffset + meshlet.triangleCount * 3 <=
                   clusters.meshletTriangles.size());

        const unsigned int *meshletVertices =
            clusters.meshletVertices.data() + meshlet.vertexOffset;
        const unsigned char *triangleData =
            clusters.meshletTriangles.data() + meshlet.triangleOffset;

        for (uint32_t triangleIndex = 0; triangleIndex < meshlet.triangleCount; triangleIndex++)
        {
            const unsigned char localIndex0 = triangleData[3 * triangleIndex];
            const unsigned char localIndex1 = triangleData[3 * triangleIndex + 1];
            const unsigned char localIndex2 = triangleData[3 * triangleIndex + 2];
            YBI_ASSERT(localIndex0 < meshlet.vertexCount);
            YBI_ASSERT(localIndex1 < meshlet.vertexCount);
            YBI_ASSERT(localIndex2 < meshlet.vertexCount);

            const TriangleKey key = {{static_cast<int>(meshletVertices[localIndex0]),
                                      static_cast<int>(meshletVertices[localIndex1]),
                                      static_cast<int>(meshletVertices[localIndex2])}};
            auto it = triangleMap.find(key);
            YBI_ASSERT(it != triangleMap.end());

            TriangleOccurrences &occurrences = it->second;
            YBI_ASSERT(occurrences.next < occurrences.sourceTriangles.size());
            triangles[outputTriangleIndex++] = occurrences.sourceTriangles[occurrences.next++];
        }
    }

    YBI_ASSERT(outputTriangleIndex == triangles.size());
}

struct TriangleNeighborRange
{
    uint32_t offset;
    uint32_t count;
};

// NOTE: neighborRanges[triangle] specifies a span in neighbors of triangles that share an edge
// with triangle
struct TriangleAdjacency
{
    View<TriangleNeighborRange> neighborRanges;
    View<uint32_t> neighbors;
};

static uint16_t CalculateEdgeKey(uint8_t i0, uint8_t i1)
{
    uint8_t lo = Min(i0, i1);
    uint8_t hi = Max(i0, i1);
    uint16_t key = lo | (static_cast<uint16_t>(hi) << 8u);
    return key;
}

static TriangleAdjacency CreateTriangleAdjacency(Arena *arena,
                                                 const MeshletClustersResult &clusters,
                                                 const Meshlet &meshlet)
{
    View<TriangleNeighborRange> triangleNeighborRanges =
        arena->PushArray<TriangleNeighborRange>(meshlet.triangleCount);

    std::unordered_multimap<uint16_t, unsigned char> edgeToClusterMap;
    const unsigned char *triangles = clusters.meshletTriangles.data() + meshlet.triangleOffset;

    // Construct edge key to triangle index map
    for (uint32_t triangleIndex = 0; triangleIndex < meshlet.triangleCount; triangleIndex++)
    {
        for (uint32_t edgeIndex = 0; edgeIndex < 3; edgeIndex++)
        {
            uint8_t index0 = static_cast<uint8_t>(triangles[3 * triangleIndex + edgeIndex]);
            uint8_t index1 =
                static_cast<uint8_t>(triangles[3 * triangleIndex + ((edgeIndex + 1) % 3)]);

            uint16_t key = CalculateEdgeKey(index0, index1);
            edgeToClusterMap.emplace(key, triangleIndex);
        }
    }

    // Calculate per triangle adjacency offsets
    uint32_t totalCount = 0;
    for (uint32_t triangleIndex = 0; triangleIndex < meshlet.triangleCount; triangleIndex++)
    {
        uint32_t adjacencyCount = 0;
        for (uint32_t edgeIndex = 0; edgeIndex < 3; edgeIndex++)
        {
            uint8_t index0 = static_cast<uint8_t>(triangles[3 * triangleIndex + edgeIndex]);
            uint8_t index1 =
                static_cast<uint8_t>(triangles[3 * triangleIndex + ((edgeIndex + 1) % 3)]);

            uint16_t key = CalculateEdgeKey(index0, index1);

            const uint32_t count = edgeToClusterMap.count(key);
            YBI_ASSERT(count > 0);
            adjacencyCount += count - 1;
        }

        TriangleNeighborRange &adjacency = triangleNeighborRanges[triangleIndex];
        adjacency.offset = totalCount;
        adjacency.count = adjacencyCount;

        totalCount += adjacencyCount;
    }

    // Fill in adjacency data
    View<uint32_t> triangleNeighbors = arena->PushArray<uint32_t>(totalCount);
    for (uint32_t triangleIndex = 0; triangleIndex < meshlet.triangleCount; triangleIndex++)
    {
        TriangleNeighborRange &adjacency = triangleNeighborRanges[triangleIndex];
        uint32_t offset = adjacency.offset;
        for (uint32_t edgeIndex = 0; edgeIndex < 3; edgeIndex++)
        {
            uint8_t index0 = static_cast<uint8_t>(triangles[3 * triangleIndex + edgeIndex]);
            uint8_t index1 =
                static_cast<uint8_t>(triangles[3 * triangleIndex + ((edgeIndex + 1) % 3)]);

            uint16_t key = CalculateEdgeKey(index0, index1);

            auto range = edgeToClusterMap.equal_range(key);
            for (auto it = range.first; it != range.second; ++it)
            {
                const uint32_t otherTriangleIndex = it->second;
                if (otherTriangleIndex != triangleIndex)
                {
                    // TODO: make sure input topology is not degenerate (e.g triangle neighboring
                    // another triangle more than once)
                    triangleNeighbors[offset++] = otherTriangleIndex;
                }
            }
        }
    }

    TriangleAdjacency adjacency;
    adjacency.neighborRanges = triangleNeighborRanges;
    adjacency.neighbors = triangleNeighbors;
    return adjacency;
}

enum TriangleStripType : int
{
    Strip_Restart,
    Strip_Edge1,
    Strip_Edge2,
    Strip_Backtrack,
};

// The (up to two) unconsumed edges of a triangle encoded in a strip (either the last or second to
// last triangle)
typedef uint32_t Triangle[3];

static void
SetRotatedTriangle(Triangle triangle, const unsigned char *triangles, uint32_t triIndex, int rot)
{
    static const int rotationsTable[] = {0, 1, 2, 0, 1};
    const int *order = rotationsTable + rot;
    triangle[0] = triangles[3 * triIndex + order[0]];
    triangle[1] = triangles[3 * triIndex + order[1]];
    triangle[2] = triangles[3 * triIndex + order[2]];
}

static int GetConnectingEdge(Triangle triangle, uint32_t a, uint32_t b, uint32_t c)
{
    const int table[] = {0, 1, 2, 0};

    // NOTE: if winding doesn't align (e.g non manifold edge), this will return -1, triggering a
    // restart
    for (uint32_t i = 0; i < 3; i++)
    {
        uint32_t e0 = triangle[table[i]];
        uint32_t e1 = triangle[table[i + 1]];

        if (a == e1 && b == e0)
        {
            return (i << 2) | 0;
        }
        else if (b == e1 && c == e0)
        {
            return (i << 2) | 1;
        }
        else if (c == e1 && a == e0)
        {
            return (i << 2) | 2;
        }
    }
    return -1;
}

static uint32_t FindMinValenceTriangle(const BitVector &triangleBitVector,
                                       const TriangleAdjacency &adjacency)
{
    uint32_t minValence = ~0u;
    uint32_t selectedTriangle = ~0u;
    for (uint32_t triangleIndex = 0; triangleIndex < adjacency.neighborRanges.size();
         triangleIndex++)
    {
        if (triangleBitVector.GetBit(triangleIndex))
            continue;

        const uint32_t count = adjacency.neighborRanges[triangleIndex].count;
        if (count < minValence)
        {
            selectedTriangle = triangleIndex;
            minValence = count;
        }
    }

    YBI_ASSERT(selectedTriangle != ~0u);
    return selectedTriangle;
}

static uint32_t RemoveTriangleAndFindMinValenceNeighbor(TriangleAdjacency &adjacency,
                                                        uint32_t triangle)
{
    uint32_t minValence = ~0u;
    uint32_t selectedTriangle = ~0u;
    uint32_t minI = ~0u;

    TriangleNeighborRange &range = adjacency.neighborRanges[triangle];
    YBI_ASSERT(range.count);
    for (uint32_t i = range.offset; i < range.offset + range.count; i++)
    {
        uint32_t neighbor = adjacency.neighbors[i];

        TriangleNeighborRange &neighborRange = adjacency.neighborRanges[neighbor];
        uint32_t neighborCount = neighborRange.count;

        bool found = false;
        for (int j = 0; j < neighborCount; j++)
        {
            if (adjacency.neighbors[neighborRange.offset + j] == triangle)
            {
                Swap(adjacency.neighbors[neighborRange.offset + j],
                     adjacency.neighbors[neighborRange.offset + neighborCount - 1]);
                neighborRange.count--;
                found = true;
                break;
            }
        }
        YBI_ASSERT(found);
        if (neighborRange.count < minValence)
        {
            selectedTriangle = neighbor;
            minValence = neighborRange.count;
            minI = i;
        }
    }
    YBI_ASSERT(selectedTriangle != ~0u);
    Swap(adjacency.neighbors[minI], adjacency.neighbors[range.offset + range.count - 1]);
    range.count--;

    return selectedTriangle;
}

static uint32_t FindOppositeEdgeNeighbor(TriangleAdjacency &adjacency,
                                         uint32_t triangle,
                                         const unsigned char *triangles,
                                         Triangle prevTriangle,
                                         int edge)
{
    YBI_ASSERT(edge == 1 || edge == 2);

    TriangleNeighborRange &range = adjacency.neighborRanges[triangle];

    uint32_t result = ~0u;
    for (uint32_t i = range.offset; i < range.offset + range.count; i++)
    {
        uint32_t neighbor = adjacency.neighbors[i];
        int code = GetConnectingEdge(prevTriangle,
                                     triangles[3 * neighbor + 0],
                                     triangles[3 * neighbor + 1],
                                     triangles[3 * neighbor + 2]);
        if (code != -1 && (code >> 2) + edge == 3)
        {
            result = neighbor;
            Swap(adjacency.neighbors[i], adjacency.neighbors[range.offset + range.count - 1]);
            range.count--;
            break;
        }
    }
    return result;
}

struct TriangleOrderingResult
{
    Array<uint32_t> triangleOrder;
    Array<int> rotations;
    Array<TriangleStripType> stripTypes;
};

static void CalculateTriangleStripOrder(TriangleAdjacency &adjacency,
                                        const MeshletClustersResult &clusters,
                                        const Meshlet &meshlet,
                                        TriangleOrderingResult &result,
                                        uint32_t clusterIndex)
{
    ScratchArena scratch;

    Array<uint32_t> triangleOrder;
    Array<int> rotations;
    Array<TriangleStripType> stripTypes;

    BitVector visitedTriangles(scratch.arena, meshlet.triangleCount);

    const unsigned char *triangles = clusters.meshletTriangles.data() + meshlet.triangleOffset;

    Triangle triangleHistory[2];
    // Current history slot represents the last triangle. The other slot (!currentHistorySlot) is
    // the second to last triangle.
    uint32_t currentHistorySlot = 0;

    for (uint32_t outputTriangleIndex = 0; outputTriangleIndex < meshlet.triangleCount;
         outputTriangleIndex++)
    {
        bool restart = true;
        // Edge 1 or 2
        if (triangleOrder.size() && adjacency.neighborRanges[triangleOrder.Last()].count)
        {
            uint32_t triangle = triangleOrder.Last();
            uint32_t nextTriangle = RemoveTriangleAndFindMinValenceNeighbor(adjacency, triangle);
            int connection = GetConnectingEdge(triangleHistory[currentHistorySlot],
                                               triangles[3 * nextTriangle + 0],
                                               triangles[3 * nextTriangle + 1],
                                               triangles[3 * nextTriangle + 2]);

            if (connection != -1)
            {
                int rotation = connection & 3;

                int connectingEdge = connection >> 2;

                // If the connecting edge is 0, there are two possible cases:
                // 1. The previous triangle was a restart. In this case, rotate the previous
                // triangle so that nextTriangle connects to edge1.
                // 2. The edge is non manifold. In this case, restart the strip.
                if (connectingEdge == 0 && stripTypes.Last() == Strip_Restart)
                {
                    // 0, 1, 2 -> 2, 0, 1
                    rotations.Last() = 2;
                    connectingEdge = 1;
                    SetRotatedTriangle(triangleHistory[currentHistorySlot],
                                       triangles,
                                       triangleOrder.Last(),
                                       rotations.Last());
                }
                if (!(connectingEdge == 0 && stripTypes.Last() != Strip_Restart))
                {
                    stripTypes.EmplaceBack(static_cast<TriangleStripType>(connectingEdge));
                    triangleOrder.EmplaceBack(nextTriangle);
                    rotations.EmplaceBack(rotation);
                    restart = false;
                }
            }
        }
        // Backtrack
        else if (triangleOrder.size() > 1 && stripTypes.Last() != Strip_Restart &&
                 stripTypes.Last() != Strip_Backtrack)
        {
            uint32_t prevTriangle = triangleOrder[triangleOrder.size() - 2];
            uint32_t nextTriangle = FindOppositeEdgeNeighbor(adjacency,
                                                             prevTriangle,
                                                             triangles,
                                                             triangleHistory[!currentHistorySlot],
                                                             stripTypes.Last());
            if (nextTriangle != ~0u)
            {
                YBI_ASSERT(visitedTriangles.GetBit(nextTriangle) == 0);
                int connection = GetConnectingEdge(triangleHistory[!currentHistorySlot],
                                                   triangles[3 * nextTriangle + 0],
                                                   triangles[3 * nextTriangle + 1],
                                                   triangles[3 * nextTriangle + 2]);
                if (connection != -1)
                {
                    int rotation = connection & 3;
                    int connectingEdge = connection >> 2;
                    YBI_ASSERT((connectingEdge == 1 || connectingEdge == 2) &&
                               connectingEdge != stripTypes.Last());

                    stripTypes.EmplaceBack(Strip_Backtrack);
                    triangleOrder.EmplaceBack(nextTriangle);
                    rotations.EmplaceBack(rotation);
                    restart = false;
                }
            }
        }

        // Restart
        if (restart)
        {
            uint32_t nextTriangle = FindMinValenceTriangle(visitedTriangles, adjacency);

            stripTypes.EmplaceBack(Strip_Restart);
            triangleOrder.EmplaceBack(nextTriangle);
            rotations.EmplaceBack(0);
        }

        visitedTriangles.SetBit(triangleOrder.Last());

        // Flips b/t 0 and 1
        currentHistorySlot = !currentHistorySlot;

        SetRotatedTriangle(triangleHistory[currentHistorySlot],
                           triangles,
                           triangleOrder.Last(),
                           rotations.Last());
    }

    result.triangleOrder = std::move(triangleOrder);
    result.rotations = std::move(rotations);
    result.stripTypes = std::move(stripTypes);
}

static std::atomic<int64_t> amountSaved{0};
static std::atomic<uint64_t> positionBytesSaved{0};
static std::atomic<uint64_t> clusterVertexCount{0};

static uint32_t CalculateBitWidth(uint32_t diff)
{
    uint32_t result = diff == 0 ? 0u : FindMostSignificantBit(diff) + 1u;
    return result;
}

struct PackedStreamInfo
{
    uint32_t base;
    uint32_t bitWidth;
    size_t byteOffset;
    size_t byteCount;
};

static size_t AppendOutputBytes(Array<uint8_t> &outputBytes, size_t byteCount)
{
    const size_t byteOffset = outputBytes.size();
    outputBytes.Resize(byteOffset + byteCount);
    if (byteCount)
    {
        memset(outputBytes.data() + byteOffset, 0, byteCount);
    }

    return byteOffset;
}

template <typename Range, typename GetValue>
static PackedStreamInfo
AppendPackedDeltas(Array<uint8_t> &outputBytes, const Range &items, GetValue getValue)
{
    YBI_ASSERT(items.size() > 0);

    uint32_t minValue = ~0u;
    uint32_t maxValue = 0;

    for (const auto &item : items)
    {
        const uint32_t value = getValue(item);
        minValue = Min(minValue, value);
        maxValue = Max(maxValue, value);
    }

    const uint32_t bitWidth = CalculateBitWidth(maxValue - minValue);
    const size_t totalBits = size_t(bitWidth) * items.size();
    const size_t byteCount = (totalBits + 7) / 8;

    const size_t byteOffset = AppendOutputBytes(outputBytes, byteCount);
    YBI_ASSERT(byteOffset <= UINT32_MAX / 8u);
    uint32_t bitOffset = static_cast<uint32_t>(byteOffset * 8u);
    for (const auto &item : items)
    {
        WriteBits(outputBytes.data(), bitOffset, getValue(item) - minValue, bitWidth);
    }

    return {minValue, bitWidth, byteOffset, byteCount};
}

// i feel like this is the biggest problem with templates. like in order for this to work with
// static

static uint32_t WriteMaterialIDs(const DataContainer *container,
                                 const Mesh &mesh,
                                 const TriangleOrderingResult &orderingInfo,
                                 const Meshlet &meshlet,
                                 const Array<int> &clusterTriangleToMeshTriangle,
                                 Array<uint8_t> &outputBytes)
{
    const BufferHandle faceMaterialHandle = mesh.GetFaceMaterialIndicesHandle();
    if (!faceMaterialHandle.IsValid())
    {
        const uint32_t materialIndex = mesh.GetMaterialIndex();
        if (materialIndex == DEVICE_CLUSTER_INVALID_MATERIAL_ID)
        {
            return DEVICE_CLUSTER_INVALID_MATERIAL_ID;
        }

        YBI_ASSERT(materialIndex < DEVICE_CLUSTER_ENCODED_INVALID_MATERIAL_ID);
        return materialIndex | 0x80000000u;
    }

    YBI_ASSERT(faceMaterialHandle.index < container->meshFaceMaterialIndices.size());
    YBI_ASSERT(clusterTriangleToMeshTriangle.size());
    const Array<int> &faceMaterialIndices =
        container->meshFaceMaterialIndices[faceMaterialHandle.index];

    Array<uint32_t> materialIndices;
    materialIndices.Reserve(meshlet.triangleCount);
    for (uint32_t outputTriangleIndex = 0; outputTriangleIndex < meshlet.triangleCount;
         outputTriangleIndex++)
    {
        const uint32_t localTriangle = orderingInfo.triangleOrder[outputTriangleIndex];
        const int meshTriangle =
            clusterTriangleToMeshTriangle[meshlet.triangleOffset / 3 + localTriangle];
        YBI_ASSERT(meshTriangle >= 0 &&
                   static_cast<size_t>(meshTriangle) < faceMaterialIndices.size());

        const int sourceMaterialIndex = faceMaterialIndices[meshTriangle];
        YBI_ASSERT(sourceMaterialIndex < 0 || static_cast<uint32_t>(sourceMaterialIndex) <
                                                  DEVICE_CLUSTER_ENCODED_INVALID_MATERIAL_ID);
        const uint32_t materialIndex = sourceMaterialIndex < 0
                                           ? DEVICE_CLUSTER_ENCODED_INVALID_MATERIAL_ID
                                           : static_cast<uint32_t>(sourceMaterialIndex);
        materialIndices.EmplaceBack(materialIndex);
    }

    YBI_ASSERT(materialIndices.size());
    const uint32_t firstMaterialIndex = materialIndices[0];
    YBI_ASSERT(firstMaterialIndex <= DEVICE_CLUSTER_ENCODED_INVALID_MATERIAL_ID);

    bool constantMaterialID = true;
    for (uint32_t materialIndex : materialIndices)
    {
        YBI_ASSERT(materialIndex <= DEVICE_CLUSTER_ENCODED_INVALID_MATERIAL_ID);
        if (materialIndex != firstMaterialIndex)
        {
            constantMaterialID = false;
        }
    }

    if (constantMaterialID)
    {
        return firstMaterialIndex | 0x80000000u;
    }

    Array<uint32_t> uniqueMaterialIDs;
    uniqueMaterialIDs.Reserve(materialIndices.size());
    Array<uint32_t> entryIndices(materialIndices.size());

    uniqueMaterialIDs.EmplaceBack(firstMaterialIndex);
    entryIndices[0] = 0;

    uint32_t differingBits = 0;
    for (size_t triangleIndex = 1; triangleIndex < materialIndices.size(); triangleIndex++)
    {
        const uint32_t materialID = materialIndices[triangleIndex];
        differingBits |= firstMaterialIndex ^ materialID;

        size_t entryIndex = 0;
        while (entryIndex < uniqueMaterialIDs.size() &&
               uniqueMaterialIDs[entryIndex] != materialID)
        {
            entryIndex++;
        }
        if (entryIndex == uniqueMaterialIDs.size())
        {
            uniqueMaterialIDs.EmplaceBack(materialID);
        }
        entryIndices[triangleIndex] = static_cast<uint32_t>(entryIndex);
    }

    YBI_ASSERT(differingBits != 0);
    const uint32_t commonMSBs = 31u - static_cast<uint32_t>(FindMostSignificantBit(differingBits));
    const uint32_t lsbBitCount = 32u - commonMSBs;
    const uint32_t commonPrefix = firstMaterialIndex >> lsbBitCount;
    const uint32_t lsbMask = (1u << lsbBitCount) - 1u;
    const uint32_t entryBitLength =
        CalculateBitWidth(static_cast<uint32_t>(uniqueMaterialIDs.size() - 1));

    YBI_ASSERT(commonMSBs > 0);
    YBI_ASSERT(uniqueMaterialIDs.size() > 1 && uniqueMaterialIDs.size() <= 31);

    const size_t materialTableNumBits = commonMSBs + lsbBitCount * uniqueMaterialIDs.size() +
                                        size_t(entryBitLength) * materialIndices.size();
    const size_t materialTableSize = (materialTableNumBits + 7u) >> 3u;
    const size_t materialTableOffset = AppendOutputBytes(outputBytes, materialTableSize);
    YBI_ASSERT(materialTableOffset < (1ull << 22u));
    YBI_ASSERT(materialTableOffset <= UINT32_MAX / 8u);

    uint32_t tableBitOffset = static_cast<uint32_t>(materialTableOffset * 8u);
    WriteBits(outputBytes.data(), tableBitOffset, commonPrefix, commonMSBs);
    for (uint32_t materialID : uniqueMaterialIDs)
    {
        WriteBits(outputBytes.data(), tableBitOffset, materialID & lsbMask, lsbBitCount);
    }
    for (uint32_t entryIndex : entryIndices)
    {
        WriteBits(outputBytes.data(), tableBitOffset, entryIndex, entryBitLength);
    }
    YBI_ASSERT(tableBitOffset == materialTableOffset * 8u + materialTableNumBits);

    return (commonMSBs - 1u) | (static_cast<uint32_t>(uniqueMaterialIDs.size()) << 5u) |
           (static_cast<uint32_t>(materialTableOffset) << 10u);
}

template <typename T>
struct PrimvarPacking;

template <>
struct PrimvarPacking<Vec2f>
{
    using PackedType = Vec2u;
    static constexpr size_t componentCount = 2;
};

template <>
struct PrimvarPacking<Vec3f>
{
    using PackedType = Vec3u;
    static constexpr size_t componentCount = 3;
};

template <>
struct PrimvarPacking<Vec4f>
{
    using PackedType = Vec4u;
    static constexpr size_t componentCount = 4;
};

template <typename T>
struct PrimvarHash
{
    using PackedType = typename PrimvarPacking<T>::PackedType;

    size_t operator()(const PackedType &primvar) const
    {
        size_t hash = 0;
        for (size_t component = 0; component < PrimvarPacking<T>::componentCount; ++component)
        {
            hash ^= std::hash<uint32_t>{}(primvar[component]) + 0x9e3779b9u + (hash << 6u) +
                    (hash >> 2u);
        }
        return hash;
    }
};

template <typename T>
static DeviceClusterAttribute WriteClusteredPrimvars(const DataContainer *container,
                                                     const Primvar &primvar,
                                                     AttributeType attributeType,
                                                     const Array<int> *primvarIndices,
                                                     const Array<unsigned int> &faceVaryingIndices,
                                                     const Array<unsigned int> &vertexIndices,
                                                     Array<uint8_t> &outputBytes)
{
    using PackedType = typename PrimvarPacking<T>::PackedType;
    constexpr size_t componentCount = PrimvarPacking<T>::componentCount;

    const std::vector<Array<T>> &arrays = container->attributes.GetArrays<T>();
    YBI_ASSERT(primvar.index < arrays.size());
    const Array<T> &primvars = arrays[primvar.index];
    const Array<unsigned int> &chosenIndices =
        attributeType == AttributeType::FaceVarying ? faceVaryingIndices : vertexIndices;
    YBI_ASSERT(chosenIndices.size());

    const Array<int> *indices = primvarIndices;

    Array<PackedType> logicalPrimvars;
    Array<PackedType> uniquePrimvars;
    Array<uint32_t> localPrimvarIndices;
    logicalPrimvars.Reserve(chosenIndices.size());
    uniquePrimvars.Reserve(chosenIndices.size());
    localPrimvarIndices.Reserve(chosenIndices.size());
    std::unordered_map<PackedType, uint32_t, PrimvarHash<T>> uniquePrimvarMap;
    uniquePrimvarMap.reserve(chosenIndices.size());

    for (unsigned int index : chosenIndices)
    {
        int primvarIndex = static_cast<int>(index);
        if (indices)
        {
            YBI_ERROR(
                index < indices->size(), "%u %zu %u\n", index, indices->size(), attributeType);
            primvarIndex = (*indices)[index];
        }
        YBI_ASSERT(primvarIndex >= 0 && static_cast<size_t>(primvarIndex) < primvars.size());

        const T &primvarValue = primvars[static_cast<size_t>(primvarIndex)];
        PackedType packedPrimvar(0u);
        for (size_t component = 0; component < componentCount; ++component)
        {
            packedPrimvar[component] = FloatAsUint(primvarValue[component]);
        }
        logicalPrimvars.EmplaceBack(packedPrimvar);

        const auto result =
            uniquePrimvarMap.emplace(packedPrimvar, static_cast<uint32_t>(uniquePrimvars.size()));
        if (result.second)
        {
            uniquePrimvars.EmplaceBack(packedPrimvar);
        }
        localPrimvarIndices.EmplaceBack(result.first->second);
    }

    YBI_ASSERT(logicalPrimvars.size() == chosenIndices.size());
    YBI_ASSERT(uniquePrimvars.size());

    PackedType minPrimvar(~0u);
    PackedType maxPrimvar(0u);
    PackedType bitWidths(0u);

    for (const PackedType &packedPrimvar : uniquePrimvars)
    {
        minPrimvar = Min(minPrimvar, packedPrimvar);
        maxPrimvar = Max(maxPrimvar, packedPrimvar);
    }

    size_t bitsPerPrimvar = 0;
    for (size_t component = 0; component < componentCount; ++component)
    {
        bitWidths[component] = CalculateBitWidth(maxPrimvar[component] - minPrimvar[component]);
        bitsPerPrimvar += bitWidths[component];
    }
    const uint32_t primvarIndexBitWidth =
        CalculateBitWidth(static_cast<uint32_t>(uniquePrimvars.size() - 1u));

    DeviceClusterAttribute attribute = {};
    for (size_t component = 0; component < componentCount; ++component)
    {
        YBI_ASSERT(bitWidths[component] <= DEVICE_CLUSTER_ATTRIBUTE_COMPONENT_BIT_WIDTH_MASK);
        attribute.bitWidths[component] = bitWidths[component];
    }

    const size_t primvarDataByteOffset =
        AppendOutputBytes(outputBytes, componentCount * sizeof(uint32_t));
    YBI_ASSERT(primvarDataByteOffset <= DEVICE_CLUSTER_ATTRIBUTE_DATA_BYTE_OFFSET_MASK);
    YBI_ASSERT(primvarDataByteOffset <= UINT32_MAX / 8u);
    attribute.dataByteOffset = static_cast<uint32_t>(primvarDataByteOffset);

    uint32_t primvarBaseBitOffset = static_cast<uint32_t>(primvarDataByteOffset * 8u);
    for (size_t component = 0; component < componentCount; ++component)
    {
        WriteBits(outputBytes.data(), primvarBaseBitOffset, minPrimvar[component], 32u);
    }

    const size_t flatPrimvarByteCount = (logicalPrimvars.size() * bitsPerPrimvar + 7u) >> 3u;
    const size_t indexedPrimvarByteCount = (uniquePrimvars.size() * bitsPerPrimvar + 7u) >> 3u;
    const size_t primvarIndexByteCount =
        (localPrimvarIndices.size() * size_t(primvarIndexBitWidth) + 7u) >> 3u;
    const bool useIndexedPrimvars =
        indexedPrimvarByteCount + primvarIndexByteCount < flatPrimvarByteCount;

    const Array<PackedType> &storedPrimvars =
        useIndexedPrimvars ? uniquePrimvars : logicalPrimvars;
    const size_t primvarByteCount = (storedPrimvars.size() * bitsPerPrimvar + 7u) >> 3u;
    const size_t primvarByteOffset = AppendOutputBytes(outputBytes, primvarByteCount);
    YBI_ASSERT(primvarByteOffset <= UINT32_MAX / 8u);
    uint32_t primvarBitOffset = static_cast<uint32_t>(primvarByteOffset * 8u);
    for (const PackedType &packedPrimvar : storedPrimvars)
    {
        for (size_t component = 0; component < componentCount; ++component)
        {
            WriteBits(outputBytes.data(),
                      primvarBitOffset,
                      packedPrimvar[component] - minPrimvar[component],
                      bitWidths[component]);
        }
    }

    if (useIndexedPrimvars)
    {
        YBI_ASSERT(primvarIndexBitWidth <= DEVICE_CLUSTER_ATTRIBUTE_INDEX_BIT_WIDTH_MASK);
        const size_t primvarIndicesByteOffset =
            AppendOutputBytes(outputBytes, primvarIndexByteCount);
        YBI_ASSERT(primvarIndicesByteOffset <= UINT32_MAX / 8u);
        YBI_ASSERT(primvarIndicesByteOffset >= primvarDataByteOffset);
        const size_t relativePrimvarIndexByteOffset =
            primvarIndicesByteOffset - primvarDataByteOffset;
        YBI_ASSERT(relativePrimvarIndexByteOffset <=
                   DEVICE_CLUSTER_ATTRIBUTE_INDEX_BYTE_OFFSET_MASK);
        attribute.indexed = true;
        attribute.indexBitWidth = primvarIndexBitWidth;
        attribute.indexByteOffset = static_cast<uint32_t>(relativePrimvarIndexByteOffset);

        uint32_t primvarIndexBitOffset = static_cast<uint32_t>(primvarIndicesByteOffset * 8u);
        for (uint32_t localPrimvarIndex : localPrimvarIndices)
        {
            WriteBits(outputBytes.data(),
                      primvarIndexBitOffset,
                      localPrimvarIndex,
                      primvarIndexBitWidth);
        }
    }

    return attribute;
}

static uint64_t PackClusterAttribute(const DeviceClusterAttribute &attribute)
{
    YBI_ASSERT(attribute.dataByteOffset <= DEVICE_CLUSTER_ATTRIBUTE_DATA_BYTE_OFFSET_MASK);

    uint64_t bits = attribute.dataByteOffset;
    for (uint32_t component = 0; component < DEVICE_CLUSTER_ATTRIBUTE_MAX_COMPONENT_COUNT;
         ++component)
    {
        YBI_ASSERT(attribute.bitWidths[component] <=
                   DEVICE_CLUSTER_ATTRIBUTE_COMPONENT_BIT_WIDTH_MASK);
        const uint32_t shift = DEVICE_CLUSTER_ATTRIBUTE_COMPONENT_BIT_WIDTH_SHIFT +
                               component * DEVICE_CLUSTER_ATTRIBUTE_COMPONENT_BIT_WIDTH_NUM_BITS;
        bits |= static_cast<uint64_t>(attribute.bitWidths[component]) << shift;
    }

    YBI_ASSERT(attribute.indexBitWidth <= DEVICE_CLUSTER_ATTRIBUTE_INDEX_BIT_WIDTH_MASK);
    YBI_ASSERT(attribute.indexByteOffset <= DEVICE_CLUSTER_ATTRIBUTE_INDEX_BYTE_OFFSET_MASK);
    bits |= static_cast<uint64_t>(attribute.indexBitWidth)
            << DEVICE_CLUSTER_ATTRIBUTE_INDEX_BIT_WIDTH_SHIFT;
    bits |= static_cast<uint64_t>(attribute.indexByteOffset)
            << DEVICE_CLUSTER_ATTRIBUTE_INDEX_BYTE_OFFSET_SHIFT;
    bits |= static_cast<uint64_t>(attribute.indexed)
            << DEVICE_CLUSTER_ATTRIBUTE_INDEXED_FLAG_SHIFT;
    return bits;
}

static uint32_t WriteClusterAttributeTable(const Array<DeviceClusterAttribute> &attributes,
                                           Array<uint8_t> &outputBytes)
{
    if (attributes.size() == 0)
    {
        return 0u;
    }

    const size_t tableByteOffset =
        AppendOutputBytes(outputBytes, attributes.size() * sizeof(uint64_t));
    YBI_ASSERT(tableByteOffset <= UINT32_MAX / 8u);

    uint32_t tableBitOffset = static_cast<uint32_t>(tableByteOffset * 8u);
    for (const DeviceClusterAttribute &attribute : attributes)
    {
        const uint64_t bits = PackClusterAttribute(attribute);
        WriteBits(outputBytes.data(), tableBitOffset, static_cast<uint32_t>(bits), 32u);
        WriteBits(outputBytes.data(), tableBitOffset, static_cast<uint32_t>(bits >> 32u), 32u);
    }

    return static_cast<uint32_t>(tableByteOffset);
}

static void CreateCompressedGeometry(const DataContainer *container,
                                     const MeshletClustersResult &clusters,
                                     const Mesh &mesh,
                                     const TriangleOrderingResult &orderingInfo,
                                     const Meshlet &meshlet,
                                     const Array<int> &clusterTriangleToMeshTriangle,
                                     Array<uint8_t> &compressedBytes,
                                     DeviceCluster &cluster)
{
    ScratchArena scratch;
    Array<unsigned char> indexBuffer;
    indexBuffer.Reserve(meshlet.triangleCount);

    static const int rotationsTable[] = {0, 1, 2, 0, 1};

    const unsigned char *triangles = clusters.meshletTriangles.data() + meshlet.triangleOffset;

    for (uint32_t triangleIndex = 0; triangleIndex < orderingInfo.triangleOrder.size();
         triangleIndex++)
    {
        uint32_t triangle = orderingInfo.triangleOrder[triangleIndex];
        TriangleStripType stripType = orderingInfo.stripTypes[triangleIndex];
        if (stripType == Strip_Restart)
        {
            const int *order = rotationsTable + orderingInfo.rotations[triangleIndex];
            indexBuffer.EmplaceBack(triangles[3 * triangle + order[0]]);
            indexBuffer.EmplaceBack(triangles[3 * triangle + order[1]]);
            indexBuffer.EmplaceBack(triangles[3 * triangle + order[2]]);
        }
        else
        {
            int index = rotationsTable[orderingInfo.rotations[triangleIndex] + 2];
            indexBuffer.EmplaceBack(triangles[3 * triangle + index]);
        }
    }

    const unsigned int *vertices = clusters.meshletVertices.data() + meshlet.vertexOffset;

    BitVector firstUse(scratch.arena, indexBuffer.size());
    Array<unsigned int> vertexIndices;
    Array<uint32_t> reuseBuffer;
    BitVector visitedVertices(scratch.arena, DEVICE_CLUSTER_MAX_VERTICES);
    Array<uint32_t> oldIndexToNewIndex(meshlet.triangleCount * 3);
    reuseBuffer.Reserve(indexBuffer.size());

    uint32_t addedVertexCount = 0;
    uint32_t currentFirstUseBit = 0;
    unsigned char maxReuseIndex = 0;
    for (unsigned char index : indexBuffer)
    {
        if (visitedVertices.GetBit(index))
        {
            currentFirstUseBit++;

            unsigned char newIndex = oldIndexToNewIndex[index];
            reuseBuffer.EmplaceBack(newIndex);
            maxReuseIndex = Max(newIndex, maxReuseIndex);
        }
        else
        {
            visitedVertices.SetBit(index);

            firstUse.SetBit(currentFirstUseBit);
            currentFirstUseBit++;
            oldIndexToNewIndex[index] = addedVertexCount;
            addedVertexCount++;

            vertexIndices.EmplaceBack(vertices[index]);
        }
    }

    Array<uint8_t> outputBytes;
    unsigned int minVertexIndex = ~0u;

    clusterVertexCount.fetch_add(static_cast<uint64_t>(vertexIndices.size()),
                                 std::memory_order_relaxed);
    uint32_t bitWidthX = 0u;
    uint32_t bitWidthY = 0u;
    uint32_t bitWidthZ = 0u;
    Vec3u positionAnchor(0);

    // Write positions
    {
        uint32_t minX = ~0u;
        uint32_t maxX = 0u;

        uint32_t minY = ~0u;
        uint32_t maxY = 0u;

        uint32_t minZ = ~0u;
        uint32_t maxZ = 0u;

        Vec3f minP(FLT_MAX);
        Vec3f maxP(-FLT_MAX);

        const Array<Vec3f> &positions = container->GetMeshPositions(mesh);

        for (unsigned int vertexIndex : vertexIndices)
        {
            Vec3f pos = positions[vertexIndex];

            pos = Abs(pos);
            minP = Min(minP, pos);
            maxP = Max(maxP, pos);

            minX = Min(minX, FloatAsUint(pos.x));
            maxX = Max(maxX, FloatAsUint(pos.x));

            minY = Min(minY, FloatAsUint(pos.y));
            maxY = Max(maxY, FloatAsUint(pos.y));

            minZ = Min(minZ, FloatAsUint(pos.z));
            maxZ = Max(maxZ, FloatAsUint(pos.z));
        }

        uint32_t deltaX = maxX - minX;
        uint32_t deltaY = maxY - minY;
        uint32_t deltaZ = maxZ - minZ;
        bitWidthX = CalculateBitWidth(deltaX);
        bitWidthY = CalculateBitWidth(deltaY);
        bitWidthZ = CalculateBitWidth(deltaZ);
        positionAnchor = Vec3u(minX, minY, minZ);

        const uint32_t verticesNumBytes =
            (vertexIndices.size() * (bitWidthX + bitWidthY + bitWidthZ + 3u) + 7u) >> 3u;

        size_t positionByteOffset = AppendOutputBytes(outputBytes, verticesNumBytes);
        // TODO IMPORTANT: WriteBits is limited to only 512 MB...
        uint32_t positionBitOffset = static_cast<uint32_t>(positionByteOffset * 8u);

        for (unsigned int vertexIndex : vertexIndices)
        {
            Vec3f pos = positions[vertexIndex];
            Vec3f absPos = Abs(pos);

            uint32_t bitsX = FloatAsUint(absPos.x) - minX;
            uint32_t bitsY = FloatAsUint(absPos.y) - minY;
            uint32_t bitsZ = FloatAsUint(absPos.z) - minZ;

            bitsX |= (FloatAsUint(pos.x) >> 31u) << bitWidthX;
            bitsY |= (FloatAsUint(pos.y) >> 31u) << bitWidthY;
            bitsZ |= (FloatAsUint(pos.z) >> 31u) << bitWidthZ;

            WriteBits(outputBytes.data(), positionBitOffset, bitsX, bitWidthX + 1u);
            WriteBits(outputBytes.data(), positionBitOffset, bitsY, bitWidthY + 1u);
            WriteBits(outputBytes.data(), positionBitOffset, bitsZ, bitWidthZ + 1u);
        }
    }

#if 0
    // Write vertex indices (used to access vertex uvs)
    if (mesh.GetUVType() == AttributeType::Vertex)
    {
        unsigned int maxVertexIndex = 0u;
        for (unsigned int vertexIndex : vertexIndices)
        {
            minVertexIndex = Min(minVertexIndex, vertexIndex);
            maxVertexIndex = Max(maxVertexIndex, vertexIndex);
        }

        const uint32_t numVertexBits = CalculateBitWidth(maxVertexIndex - minVertexIndex);
        const uint32_t vertexIndicesNumBytes = ((numVertexBits * vertexIndices.size()) + 7u) >> 3u;

        const size_t byteOffset = AppendOutputBytes(outputBytes, vertexIndicesNumBytes);
        YBI_ASSERT(byteOffset <= UINT32_MAX / 8u);
        uint32_t vertexBitOffset = static_cast<uint32_t>(byteOffset * 8u);

        for (unsigned int vertexIndex : vertexIndices)
        {
            WriteBits(
                outputBytes.data(), vertexBitOffset, vertexIndex - minVertexIndex, numVertexBits);
        }
    }
#endif

    Array<unsigned int> faceVaryingIndices;
    faceVaryingIndices.Reserve(3 * meshlet.triangleCount);

    bool hasFaceVaryingAttributes = mesh.GetNormalType() == AttributeType::FaceVarying ||
                                    mesh.GetUVType() == AttributeType::FaceVarying;
    for (const Primvar &primvar : mesh.GetPrimvars())
    {
        hasFaceVaryingAttributes |= primvar.type == AttributeType::FaceVarying;
    }

    // TODO: also need this if there are face varying attributes
    if (hasFaceVaryingAttributes)
    {
        static const int rotationsTable[] = {0, 1, 2, 0, 1};
        YBI_ASSERT(clusterTriangleToMeshTriangle.size());
        for (uint32_t outputTriangleIndex = 0; outputTriangleIndex < meshlet.triangleCount;
             outputTriangleIndex++)
        {
            const uint32_t localTriangle = orderingInfo.triangleOrder[outputTriangleIndex];
            const int triangle =
                clusterTriangleToMeshTriangle[meshlet.triangleOffset / 3 + localTriangle];
            const int rotation = orderingInfo.rotations[outputTriangleIndex];
            YBI_ASSERT(rotation >= 0 && rotation < 3);
            const int *order = rotationsTable + rotation;
            faceVaryingIndices.EmplaceBack(3 * triangle + order[0]);
            faceVaryingIndices.EmplaceBack(3 * triangle + order[1]);
            faceVaryingIndices.EmplaceBack(3 * triangle + order[2]);
        }
    }

    Vec2u normalBitWidths(0u);
    Vec2u normalBase(0u);
    uint32_t normalBitOffset = 0u;

    // Write octahedral normals
    if (mesh.GetNormalType() != AttributeType::None)
    {
        const Array<Vec3f> &normals = container->GetMeshNormals(mesh);
        const AttributeType normalType = mesh.GetNormalType();
        const Array<unsigned int> &chosenIndices =
            normalType == AttributeType::FaceVarying ? faceVaryingIndices : vertexIndices;

        uint16_t minOctahedralX = 0xffff;
        uint16_t minOctahedralY = 0xffff;
        uint16_t maxOctahedralX = 0;
        uint16_t maxOctahedralY = 0;

        Array<uint32_t> octahedralNormals;
        octahedralNormals.Reserve(chosenIndices.size());
        for (unsigned int index : chosenIndices)
        {
            Vec3f normal = normals[index];
            // TODO: properly handle invalid normals
            if (LengthSquared(normal) < 1e-20f)
            {
                normal = Vec3f(0.f, 0.f, 1.f);
            }
            uint32_t octahedralNormal = EncodeOctahedral(Normalize(normal));
            octahedralNormals.EmplaceBack(octahedralNormal);

            const uint16_t octX = octahedralNormal & 0xffff;
            const uint16_t octY = octahedralNormal >> 16;

            minOctahedralX = Min(minOctahedralX, octX);
            maxOctahedralX = Max(maxOctahedralX, octX);
            minOctahedralY = Min(minOctahedralY, octY);
            maxOctahedralY = Max(maxOctahedralY, octY);
        }

        const uint32_t numOctBitsX =
            CalculateBitWidth(static_cast<uint32_t>(maxOctahedralX - minOctahedralX));
        const uint32_t numOctBitsY =
            CalculateBitWidth(static_cast<uint32_t>(maxOctahedralY - minOctahedralY));
        normalBitWidths = Vec2u(numOctBitsX, numOctBitsY);
        normalBase = Vec2u(minOctahedralX, minOctahedralY);

        const size_t normalBitCount = size_t(numOctBitsX + numOctBitsY) * octahedralNormals.size();
        const size_t normalByteCount = (normalBitCount + 7u) >> 3u;

        const size_t normalByteOffset = AppendOutputBytes(outputBytes, normalByteCount);
        YBI_ASSERT(normalByteOffset <= UINT32_MAX / 8u);
        normalBitOffset = static_cast<uint32_t>(normalByteOffset * 8u);
        uint32_t normalBitWriteOffset = normalBitOffset;
        for (uint32_t octahedralNormal : octahedralNormals)
        {
            const uint32_t octX = octahedralNormal & 0xffff;
            const uint32_t octY = octahedralNormal >> 16;

            WriteBits(
                outputBytes.data(), normalBitWriteOffset, octX - minOctahedralX, numOctBitsX);
            WriteBits(
                outputBytes.data(), normalBitWriteOffset, octY - minOctahedralY, numOctBitsY);
        }
    }

    Vec2u uvBitWidths(0u);
    Vec2u uvBase(0u);
    uint32_t uvBitOffset = 0u;
    uint32_t uvIndexBitOffset = 0u;
    uint32_t uvInfo = 0u;

    // Write flat or cluster-local indexed UVs, whichever is smaller.
    if (mesh.GetUVType() != AttributeType::None)
    {
        const BufferHandle uvHandle = mesh.GetUVsHandle();
        YBI_ASSERT(uvHandle.IsValid());
        YBI_ASSERT(uvHandle.index < container->meshUVs.size());
        const Array<Vec2f> &uvs = container->meshUVs[uvHandle.index];

        const AttributeType uvType = mesh.GetUVType();
        const Array<unsigned int> &chosenIndices =
            uvType == AttributeType::FaceVarying ? faceVaryingIndices : vertexIndices;
        YBI_ASSERT(chosenIndices.size());

        const BufferHandle uvIndicesHandle = mesh.GetUVIndicesHandle();
        const Array<int> *uvIndices =
            uvIndicesHandle.IsValid() ? &container->GetMeshUVIndices(mesh) : nullptr;

        Array<Vec2u> logicalUVs;
        Array<Vec2u> uniqueUVs;
        Array<uint32_t> localUVIndices;
        logicalUVs.Reserve(chosenIndices.size());
        uniqueUVs.Reserve(chosenIndices.size());
        localUVIndices.Reserve(chosenIndices.size());
        std::unordered_map<uint64_t, uint32_t> uniqueUVMap;
        uniqueUVMap.reserve(chosenIndices.size());

        for (unsigned int index : chosenIndices)
        {
            int uvIndex = static_cast<int>(index);
            if (uvIndices)
            {
                YBI_ASSERT(index < uvIndices->size());
                uvIndex = (*uvIndices)[index];
            }
            YBI_ASSERT(uvIndex >= 0 && static_cast<size_t>(uvIndex) < uvs.size());

            const Vec2f uv = uvs[static_cast<size_t>(uvIndex)];
            const uint32_t uvX = FloatAsUint(uv.x);
            const uint32_t uvY = FloatAsUint(uv.y);
            const Vec2u uvBits(uvX, uvY);
            logicalUVs.EmplaceBack(uvBits);

            const uint64_t key = static_cast<uint64_t>(uvX) | (static_cast<uint64_t>(uvY) << 32u);
            const auto result = uniqueUVMap.emplace(key, static_cast<uint32_t>(uniqueUVs.size()));
            if (result.second)
            {
                uniqueUVs.EmplaceBack(uvBits);
            }
            localUVIndices.EmplaceBack(result.first->second);
        }

        YBI_ASSERT(logicalUVs.size() == chosenIndices.size());
        YBI_ASSERT(uniqueUVs.size());

        uint32_t minUvX = ~0u;
        uint32_t minUvY = ~0u;
        uint32_t maxUvX = 0u;
        uint32_t maxUvY = 0u;
        for (const Vec2u &uv : uniqueUVs)
        {
            minUvX = Min(minUvX, uv.x);
            maxUvX = Max(maxUvX, uv.x);
            minUvY = Min(minUvY, uv.y);
            maxUvY = Max(maxUvY, uv.y);
        }

        uvBitWidths =
            Vec2u(CalculateBitWidth(maxUvX - minUvX), CalculateBitWidth(maxUvY - minUvY));
        uvBase = Vec2u(minUvX, minUvY);
        const size_t bitsPerUV = size_t(uvBitWidths.x) + uvBitWidths.y;
        const uint32_t uvIndexBitWidth =
            CalculateBitWidth(static_cast<uint32_t>(uniqueUVs.size() - 1u));

        const size_t flatUVByteCount = (logicalUVs.size() * bitsPerUV + 7u) >> 3u;
        const size_t indexedUVByteCount = (uniqueUVs.size() * bitsPerUV + 7u) >> 3u;
        const size_t uvIndexByteCount =
            (localUVIndices.size() * size_t(uvIndexBitWidth) + 7u) >> 3u;
        const bool useIndexedUVs = indexedUVByteCount + uvIndexByteCount < flatUVByteCount;

        const Array<Vec2u> &storedUVs = useIndexedUVs ? uniqueUVs : logicalUVs;
        const size_t uvByteCount = (storedUVs.size() * bitsPerUV + 7u) >> 3u;
        const size_t uvByteOffset = AppendOutputBytes(outputBytes, uvByteCount);
        YBI_ASSERT(uvByteOffset <= UINT32_MAX / 8u);
        uvBitOffset = static_cast<uint32_t>(uvByteOffset * 8u);
        uint32_t uvBitWriteOffset = uvBitOffset;
        for (const Vec2u &uv : storedUVs)
        {
            WriteBits(outputBytes.data(), uvBitWriteOffset, uv.x - minUvX, uvBitWidths.x);
            WriteBits(outputBytes.data(), uvBitWriteOffset, uv.y - minUvY, uvBitWidths.y);
        }

        if (useIndexedUVs)
        {
            YBI_ASSERT(uvIndexBitWidth <= DEVICE_CLUSTER_UV_INDEX_BIT_WIDTH_MASK);
            const size_t uvIndicesByteOffset = AppendOutputBytes(outputBytes, uvIndexByteCount);
            YBI_ASSERT(uvIndicesByteOffset <= UINT32_MAX / 8u);
            uvIndexBitOffset = static_cast<uint32_t>(uvIndicesByteOffset * 8u);
            uint32_t uvIndexBitWriteOffset = uvIndexBitOffset;
            for (uint32_t localUVIndex : localUVIndices)
            {
                WriteBits(
                    outputBytes.data(), uvIndexBitWriteOffset, localUVIndex, uvIndexBitWidth);
            }

            uvInfo = DEVICE_CLUSTER_UV_INDEXED_FLAG |
                     (uvIndexBitWidth << DEVICE_CLUSTER_UV_INDEX_BIT_WIDTH_SHIFT);
        }
    }

    // Write primvars
    Array<DeviceClusterAttribute> clusterAttributes;
    clusterAttributes.Reserve(mesh.GetPrimvars().size());
    for (const Primvar &primvar : mesh.GetPrimvars())
    {
        DeviceClusterAttribute clusterAttribute = {};
        if (IsPrimvarArray(primvar.valueType))
        {
            const Array<int> *primvarIndices =
                primvar.HasIndices() ? &container->attributes.GetIndices(primvar.indicesIndex)
                                     : nullptr;
            switch (primvar.valueType)
            {
                case PrimvarValueType::ArrayVec2f:
                {
                    clusterAttribute = WriteClusteredPrimvars<Vec2f>(container,
                                                                     primvar,
                                                                     primvar.type,
                                                                     primvarIndices,
                                                                     faceVaryingIndices,
                                                                     vertexIndices,
                                                                     outputBytes);
                    break;
                }
                case PrimvarValueType::ArrayVec3f:
                {
                    clusterAttribute = WriteClusteredPrimvars<Vec3f>(container,
                                                                     primvar,
                                                                     primvar.type,
                                                                     primvarIndices,
                                                                     faceVaryingIndices,
                                                                     vertexIndices,
                                                                     outputBytes);
                    break;
                }
                case PrimvarValueType::ArrayVec4f:
                {
                    clusterAttribute = WriteClusteredPrimvars<Vec4f>(container,
                                                                     primvar,
                                                                     primvar.type,
                                                                     primvarIndices,
                                                                     faceVaryingIndices,
                                                                     vertexIndices,
                                                                     outputBytes);
                    break;
                }
                default:
                {
                    YBI_ASSERT(0);
                    break; // scalar primvars are not clustered yet
                }
            }
        }
        clusterAttributes.EmplaceBack(clusterAttribute);
    }

    BitVector backtrack(scratch.arena, meshlet.triangleCount);
    BitVector restart(scratch.arena, meshlet.triangleCount);
    BitVector edge(scratch.arena, meshlet.triangleCount);

    uint32_t restartHighBitPerDword[DEVICE_CLUSTER_CONTROL_DWORD_COUNT] = {0, 0, 0, 0};
    uint32_t restartCountPerDword[DEVICE_CLUSTER_CONTROL_DWORD_COUNT] = {1, 0, 0, 0};
    int32_t edge1HighBitPerDword[DEVICE_CLUSTER_CONTROL_DWORD_COUNT] = {-1, -1, -1, -1};
    int32_t edge2HighBitPerDword[DEVICE_CLUSTER_CONTROL_DWORD_COUNT] = {-1, -1, -1, -1};
    {

        int32_t prevEdge1HighBitPerDword = 0;
        int32_t prevEdge2HighBitPerDword = 0;

        YBI_ASSERT(meshlet.triangleCount > 0 &&
                   meshlet.triangleCount <= DEVICE_CLUSTER_MAX_TRIANGLES);
        YBI_ASSERT(orderingInfo.stripTypes.size() == meshlet.triangleCount);
        YBI_ASSERT(orderingInfo.stripTypes[0] == Strip_Restart);

        restart.SetBit(0);

        for (uint32_t i = 1; i < meshlet.triangleCount; i++)
        {
            const uint32_t dwordIndex = i >> 5;

            const TriangleStripType stripType = orderingInfo.stripTypes[i];
            if (stripType == Strip_Backtrack)
            {
                YBI_ASSERT(i > 1);
                const TriangleStripType prevType = orderingInfo.stripTypes[i - 1];
                YBI_ASSERT(prevType == Strip_Edge1 || prevType == Strip_Edge2);
                backtrack.SetBit(i);
                const uint32_t prevDwordIndex = (i - 1) >> 5;
                if (prevType == Strip_Edge2)
                {
                    YBI_ASSERT(edge2HighBitPerDword[prevDwordIndex] ==
                               static_cast<int32_t>(i - 1));
                    edge1HighBitPerDword[prevDwordIndex] = static_cast<int32_t>(i - 1);
                    edge2HighBitPerDword[prevDwordIndex] = prevEdge2HighBitPerDword;
                    edge.SetBit(i);
                }
                else
                {
                    YBI_ASSERT(edge1HighBitPerDword[prevDwordIndex] ==
                               static_cast<int32_t>(i - 1));
                    edge2HighBitPerDword[prevDwordIndex] = static_cast<int32_t>(i - 1);
                    edge1HighBitPerDword[prevDwordIndex] = prevEdge1HighBitPerDword;
                }
            }
            else if (stripType == Strip_Edge1)
            {
                edge.SetBit(i);
                prevEdge1HighBitPerDword = edge1HighBitPerDword[dwordIndex];
                edge1HighBitPerDword[dwordIndex] = static_cast<int32_t>(i);
            }
            else if (stripType == Strip_Edge2)
            {
                prevEdge2HighBitPerDword = edge2HighBitPerDword[dwordIndex];
                edge2HighBitPerDword[dwordIndex] = static_cast<int32_t>(i);
            }
            else
            {
                YBI_ASSERT(stripType == Strip_Restart);
                restart.SetBit(i);
                restartHighBitPerDword[dwordIndex] = i;
                restartCountPerDword[dwordIndex]++;
            }
        }

        for (uint32_t i = 1; i < DEVICE_CLUSTER_CONTROL_DWORD_COUNT; i++)
        {
            restartCountPerDword[i] += restartCountPerDword[i - 1];
            restartHighBitPerDword[i] =
                Max(restartHighBitPerDword[i], restartHighBitPerDword[i - 1]);
            edge1HighBitPerDword[i] = Max(edge1HighBitPerDword[i], edge1HighBitPerDword[i - 1]);
            edge2HighBitPerDword[i] = Max(edge2HighBitPerDword[i], edge2HighBitPerDword[i - 1]);
        }
    }

    // Write ctrl bits, reuse, first use
    const uint32_t numIndexBits = CalculateBitWidth(static_cast<uint32_t>(maxReuseIndex));
    uint32_t controlBitOffset = 0;
    uint32_t reuseBitOffset = 0;
    uint32_t firstUseBitOffset = 0;
    {
        const uint32_t ctrlBitSize = ((meshlet.triangleCount + 31) >> 5u) * 12u;

        uint32_t bitStreamSize =
            ((numIndexBits * reuseBuffer.size() + currentFirstUseBit + 7) >> 3) + ctrlBitSize;

        const size_t initialOffset = AppendOutputBytes(outputBytes, bitStreamSize);
        YBI_ASSERT(initialOffset <= UINT32_MAX / 8u);
        uint32_t bitOffset = static_cast<uint32_t>(initialOffset * 8u);
        controlBitOffset = bitOffset;

        // Write ctrl bits
        int numRemainingTriangles = meshlet.triangleCount;
        uint32_t dwordIndex = 0;
        while (numRemainingTriangles > 0)
        {
            WriteBits(outputBytes.data(), bitOffset, restart.bits[dwordIndex], 32);
            WriteBits(outputBytes.data(), bitOffset, edge.bits[dwordIndex], 32);
            WriteBits(outputBytes.data(), bitOffset, backtrack.bits[dwordIndex], 32);
            dwordIndex++;
            numRemainingTriangles -= 32;
        }

        // Write reuse
        reuseBitOffset = bitOffset;
        for (uint32_t reuseIndex : reuseBuffer)
        {
            WriteBits(outputBytes.data(), bitOffset, reuseIndex, numIndexBits);
        }

        // Write first use
        firstUseBitOffset = bitOffset;
        uint32_t firstBitWriteOffset = 0;
        uint32_t firstUseBitWriteSize = Min(currentFirstUseBit - firstBitWriteOffset, 32u);
        while (firstBitWriteOffset != currentFirstUseBit)
        {
            WriteBits(outputBytes.data(),
                      bitOffset,
                      firstUse.bits[firstBitWriteOffset >> 5],
                      firstUseBitWriteSize);
            firstBitWriteOffset += firstUseBitWriteSize;
            firstUseBitWriteSize = Min(currentFirstUseBit - firstBitWriteOffset, 32u);
        }
    }

    const uint32_t packedMaterialData = WriteMaterialIDs(
        container, mesh, orderingInfo, meshlet, clusterTriangleToMeshTriangle, outputBytes);
    const uint32_t attributeTableOffset =
        WriteClusterAttributeTable(clusterAttributes, outputBytes);

    uint32_t headerOffset = 0;
    cluster = {};
    cluster.posBitWidths = Vec3u(bitWidthX + 1u, bitWidthY + 1u, bitWidthZ + 1u);
    cluster.anchor = positionAnchor;
    cluster.normalBitWidths = normalBitWidths;
    cluster.normalBase = normalBase;
    cluster.normalBitOffset = normalBitOffset;
    cluster.uvBitWidths = uvBitWidths;
    cluster.uvBase = uvBase;
    cluster.uvBitOffset = uvBitOffset;
    cluster.uvIndexBitOffset = uvIndexBitOffset;
    cluster.uvInfo = uvInfo;
    cluster.materialInfo = packedMaterialData;
    cluster.triangleCount = meshlet.triangleCount;
    cluster.indexBitWidth = numIndexBits;
    cluster.controlBitOffset = controlBitOffset;
    cluster.reuseBitOffset = reuseBitOffset;
    cluster.firstUseBitOffset = firstUseBitOffset;
    cluster.attributeTableOffset = attributeTableOffset;

    for (uint32_t i = 0; i < DEVICE_CLUSTER_CONTROL_DWORD_COUNT; i++)
    {
        cluster.restartCountPerDword[i] = restartCountPerDword[i];
        cluster.restartHighBitPerDword[i] = static_cast<int32_t>(restartHighBitPerDword[i]);
        cluster.edge1HighBitPerDword[i] = edge1HighBitPerDword[i];
        cluster.edge2HighBitPerDword[i] = edge2HighBitPerDword[i];
    }

    const uint32_t firstUseDwordCount = (currentFirstUseBit + 31u) >> 5u;
    YBI_ASSERT(firstUseDwordCount <= DEVICE_CLUSTER_FIRST_USE_DWORD_COUNT);
    uint32_t firstUseCount = 0;
    for (uint32_t i = 0; i < firstUseDwordCount; i++)
    {
        cluster.firstUseCountBeforeDword[i] = firstUseCount;
        firstUseCount += PopCount(firstUse.bits[i]);
    }
    YBI_ASSERT(firstUseCount == vertexIndices.size());

#if 0//defined(WITH_CUDA)
    ValidateCompressedClusterCUDA(outputBytes,
                                  cluster,
                                  *container,
                                  mesh,
                                  clusters,
                                  meshlet,
                                  vertexIndices,
                                  faceVaryingIndices,
                                  clusterTriangleToMeshTriangle,
                                  orderingInfo.triangleOrder,
                                  orderingInfo.rotations,
                                  oldIndexToNewIndex);
#endif

    compressedBytes = std::move(outputBytes);
}

void CompressClusters(MeshletClustersResult &clusters,
                      const Mesh &mesh,
                      DataContainer *container,
                      Array<DeviceCluster> &deviceClusters,
                      size_t clusterOffset)
{
    Array<int> triangles;

    if (mesh.GetNormalType() == AttributeType::FaceVarying ||
        mesh.GetUVType() == AttributeType::FaceVarying ||
        mesh.GetFaceMaterialIndicesHandle().IsValid())
    {
        BuildMeshletTriangleBuffer(clusters, container->GetMeshIndices(mesh), triangles);
    }

    YBI_ASSERT(clusterOffset + clusters.meshlets.size() <=
               container->compressedClusterBytes.size());
    YBI_ASSERT(deviceClusters.size() == clusters.meshlets.size());

    tbb::parallel_for(size_t(0), clusters.meshlets.size(), [&](size_t clusterIndex) {
        ScratchArena scratch;
        const Meshlet &meshlet = clusters.meshlets[clusterIndex];

        TriangleAdjacency adjacency = CreateTriangleAdjacency(scratch.arena, clusters, meshlet);

        TriangleOrderingResult orderingInfo;
        CalculateTriangleStripOrder(
            adjacency, clusters, meshlet, orderingInfo, static_cast<uint32_t>(clusterIndex));

        CreateCompressedGeometry(container,
                                 clusters,
                                 mesh,
                                 orderingInfo,
                                 meshlet,
                                 triangles,
                                 container->compressedClusterBytes[clusterOffset + clusterIndex],
                                 deviceClusters[clusterIndex]);

        // The cluster BLAS must expose the same primitive and corner order as the
        // compressed stream so that primitive IDs and barycentrics address the
        // corresponding face-varying attributes.
        static const int rotationsTable[] = {0, 1, 2, 0, 1};
        Array<unsigned char> reorderedTriangles(meshlet.triangleCount * 3u);
        const unsigned char *sourceTriangles =
            clusters.meshletTriangles.data() + meshlet.triangleOffset;
        for (uint32_t outputTriangleIndex = 0; outputTriangleIndex < meshlet.triangleCount;
             outputTriangleIndex++)
        {
            const uint32_t sourceTriangleIndex = orderingInfo.triangleOrder[outputTriangleIndex];
            const int rotation = orderingInfo.rotations[outputTriangleIndex];
            YBI_ASSERT(rotation >= 0 && rotation < 3);
            const int *order = rotationsTable + rotation;

            reorderedTriangles[3u * outputTriangleIndex + 0u] =
                sourceTriangles[3u * sourceTriangleIndex + order[0]];
            reorderedTriangles[3u * outputTriangleIndex + 1u] =
                sourceTriangles[3u * sourceTriangleIndex + order[1]];
            reorderedTriangles[3u * outputTriangleIndex + 2u] =
                sourceTriangles[3u * sourceTriangleIndex + order[2]];
        }

        unsigned char *destinationTriangles =
            clusters.meshletTriangles.data() + meshlet.triangleOffset;
        for (size_t index = 0; index < reorderedTriangles.size(); index++)
        {
            destinationTriangles[index] = reorderedTriangles[index];
        }
    });
}

void Help()
{
    printf("total position bytes saved: %llu\n",
           static_cast<unsigned long long>(positionBytesSaved.load(std::memory_order_relaxed)));
    printf("total cluster vertices: %llu\n",
           static_cast<unsigned long long>(clusterVertexCount.load(std::memory_order_relaxed)));
    printf("total bytes saved: %zu\n", amountSaved.load());
}

} // namespace ybi
