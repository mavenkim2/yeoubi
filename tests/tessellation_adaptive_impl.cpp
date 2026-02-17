#include "bvh/usd_subdiv_select.h"
#include "bvh/usd_camera_utils.h"

#include <pxr/usd/usd/stage.h>

#include <opensubdiv/far/patchMap.h>
#include <opensubdiv/far/patchTable.h>
#include <opensubdiv/far/patchTableFactory.h>
#include <opensubdiv/far/primvarRefiner.h>
#include <opensubdiv/far/ptexIndices.h>
#include <opensubdiv/far/topologyDescriptor.h>
#include <opensubdiv/sdc/types.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <limits>
#include <string>
#include <unordered_map>
#include <vector>

namespace fs = std::filesystem;
namespace Far = OpenSubdiv::Far;
namespace Sdc = OpenSubdiv::Sdc;

template <typename T>
struct OsdData
{
    T value;

    void Clear()
    {
        value = T(0.0f);
    }
    void AddWithWeight(const OsdData &src, float w)
    {
        value = value + src.value * w;
    }
};

struct TessMesh
{
    std::vector<pxr::GfVec3f> positions;
    std::vector<int> indices;
};

struct CameraPreset
{
    const char *name;
    pxr::GfVec3f origin;
    pxr::GfVec3f target;
    pxr::GfVec3f up;
    float fovYDegrees;
    int width;
    int height;
    float distanceToTarget;
};

struct EdgeKey
{
    int faceA;
    int edgeA;
    int faceB;
    int edgeB;

    bool operator==(const EdgeKey &other) const
    {
        return faceA == other.faceA && edgeA == other.edgeA && faceB == other.faceB &&
               edgeB == other.edgeB;
    }
};

struct EdgeKeyHash
{
    size_t operator()(const EdgeKey &key) const
    {
        size_t h = 1469598103934665603ull;
        auto hashMix = [&](int value) {
            h ^= size_t(uint32_t(value));
            h *= 1099511628211ull;
        };
        hashMix(key.faceA);
        hashMix(key.edgeA);
        hashMix(key.faceB);
        hashMix(key.edgeB);
        return h;
    }
};

struct EdgeVertexKey
{
    EdgeKey edgeKey;
    int segment;

    bool operator==(const EdgeVertexKey &other) const
    {
        return edgeKey == other.edgeKey && segment == other.segment;
    }
};

struct EdgeVertexKeyHash
{
    size_t operator()(const EdgeVertexKey &key) const
    {
        EdgeKeyHash edgeHash = {};
        size_t h = edgeHash(key.edgeKey);
        h ^= size_t(uint32_t(key.segment)) + 0x9e3779b9 + (h << 6) + (h >> 2);
        return h;
    }
};

struct EdgeAdjacency
{
    std::array<int, 4> adjacentFace = {{-1, -1, -1, -1}};
    std::array<int, 4> adjacentEdge = {{-1, -1, -1, -1}};
};

struct AdaptiveSettings
{
    int numEdgeSamples = 4;
    float pixelSpacing = 1.0f;
    int minRate = 1;
    int maxRate = 32;
};

static void BuildCreasePairs(const SelectedSubdivMesh &source,
                             std::vector<int> &creaseVertexPairsOut,
                             std::vector<float> &creaseSharpnessOut)
{
    creaseVertexPairsOut.clear();
    creaseSharpnessOut.clear();

    size_t indexOffset = 0;
    const size_t creaseCount = std::min(source.creaseLengths.size(), source.creaseSharpnesses.size());
    for (size_t crease = 0; crease < creaseCount; crease++)
    {
        const int length = source.creaseLengths[crease];
        if (length < 2)
        {
            indexOffset += size_t(std::max(0, length));
            continue;
        }
        if (indexOffset + size_t(length) > source.creaseIndices.size())
        {
            break;
        }

        const float sharpness = source.creaseSharpnesses[crease];
        for (int segment = 0; segment < length - 1; segment++)
        {
            const int i0 = source.creaseIndices[indexOffset + size_t(segment)];
            const int i1 = source.creaseIndices[indexOffset + size_t(segment + 1)];
            creaseVertexPairsOut.push_back(i0);
            creaseVertexPairsOut.push_back(i1);
            creaseSharpnessOut.push_back(sharpness);
        }
        indexOffset += size_t(length);
    }
}

template <typename T>
static std::vector<OsdData<T>> InterpolateVertex(const std::vector<T> &array,
                                                 const Far::TopologyRefiner *refiner,
                                                 const Far::PatchTable *patchTable)
{
    const int numLevels = refiner->GetNumLevels();
    const int numVertices = refiner->GetNumVerticesTotal();
    const int numLocalPoints = patchTable->GetNumLocalPoints();

    std::vector<OsdData<T>> values(size_t(numVertices + numLocalPoints));
    for (size_t i = 0; i < array.size(); i++)
    {
        values[i].value = array[i];
    }

    OsdData<T> *src = values.data();
    for (int level = 1; level < numLevels; level++)
    {
        OsdData<T> *dst = src + refiner->GetLevel(level - 1).GetNumVertices();
        Far::PrimvarRefiner(*refiner).Interpolate(level, src, dst);
        src = dst;
    }
    patchTable->ComputeLocalPointValues(values.data(), values.data() + numVertices);
    return values;
}

static pxr::GfVec3f EvaluatePosition(const Far::PatchMap &patchMap,
                                     const Far::PatchTable &patchTable,
                                     const std::vector<OsdData<pxr::GfVec3f>> &positions,
                                     int ptexFace,
                                     float u,
                                     float v)
{
    const Far::PatchTable::PatchHandle *handle = patchMap.FindPatch(ptexFace, u, v);
    if (!handle)
    {
        return pxr::GfVec3f(0.0f);
    }

    Far::ConstIndexArray cvIndices = patchTable.GetPatchVertices(*handle);
    std::vector<float> pWeights(size_t(cvIndices.size()));
    patchTable.EvaluateBasis(*handle, u, v, pWeights.data());

    pxr::GfVec3f p(0.0f);
    for (int cv = 0; cv < cvIndices.size(); cv++)
    {
        p += positions[size_t(cvIndices[cv])].value * pWeights[size_t(cv)];
    }
    return p;
}

static pxr::GfVec3f Normalize(const pxr::GfVec3f &v)
{
    const float len = v.GetLength();
    if (len <= 1e-8f)
    {
        return pxr::GfVec3f(0.0f, 0.0f, 1.0f);
    }
    return v / len;
}

static pxr::GfVec3f ToScreen(const CameraPreset &camera, const pxr::GfVec3f &point)
{
    const pxr::GfVec3f forward = Normalize(camera.target - camera.origin);
    const pxr::GfVec3f right = Normalize(pxr::GfCross(forward, camera.up));
    const pxr::GfVec3f up = Normalize(pxr::GfCross(right, forward));

    const pxr::GfVec3f rel = point - camera.origin;
    const float x = pxr::GfDot(rel, right);
    const float y = pxr::GfDot(rel, up);
    const float z = std::max(pxr::GfDot(rel, forward), 1e-4f);

    const float aspect = float(camera.width) / float(camera.height);
    const float tanHalfFovY = tanf(0.5f * camera.fovYDegrees * float(M_PI / 180.0));
    const float tanHalfFovX = tanHalfFovY * aspect;

    const float ndcX = x / (z * tanHalfFovX);
    const float ndcY = y / (z * tanHalfFovY);
    const float sx = (ndcX * 0.5f + 0.5f) * float(camera.width);
    const float sy = (1.0f - (ndcY * 0.5f + 0.5f)) * float(camera.height);
    return pxr::GfVec3f(sx, sy, z);
}

static pxr::GfVec2f EdgeUV(int edge, float s)
{
    switch (edge)
    {
        case 0:
            return pxr::GfVec2f(s, 0.0f);
        case 1:
            return pxr::GfVec2f(1.0f, s);
        case 2:
            return pxr::GfVec2f(1.0f - s, 1.0f);
        default:
            return pxr::GfVec2f(0.0f, 1.0f - s);
    }
}

static int ComputeEdgeRateTmax(const Far::PatchMap &patchMap,
                               const Far::PatchTable &patchTable,
                               const std::vector<OsdData<pxr::GfVec3f>> &positions,
                               int ptexFace,
                               int edge,
                               const CameraPreset &camera,
                               const AdaptiveSettings &settings)
{
    const int sampleCount = std::max(2, settings.numEdgeSamples);
    float maxLi = 0.0f;
    float sumLi = 0.0f;

    const pxr::GfVec2f uvStart = EdgeUV(edge, 0.0f);
    const pxr::GfVec3f pStart =
        EvaluatePosition(patchMap, patchTable, positions, ptexFace, uvStart[0], uvStart[1]);
    pxr::GfVec3f screenPrev = ToScreen(camera, pStart);

    for (int i = 1; i < sampleCount; i++)
    {
        const float s = float(i) / float(sampleCount - 1);
        const pxr::GfVec2f uv = EdgeUV(edge, s);
        const pxr::GfVec3f p =
            EvaluatePosition(patchMap, patchTable, positions, ptexFace, uv[0], uv[1]);
        const pxr::GfVec3f screenP = ToScreen(camera, p);
        const float dx = screenP[0] - screenPrev[0];
        const float dy = screenP[1] - screenPrev[1];
        const float li = sqrtf(dx * dx + dy * dy);
        sumLi += li;
        maxLi = std::max(maxLi, li);
        screenPrev = screenP;
    }

    const int tMin = int(ceilf(sumLi / settings.pixelSpacing));
    const int tMax = int(ceilf(float(sampleCount - 1) * maxLi / settings.pixelSpacing));
    (void)tMin;
    return std::clamp(std::max(1, tMax), settings.minRate, settings.maxRate);
}

static std::array<CameraPreset, 2> BuildCameraPresets(const SelectedSubdivMesh &selected,
                                                      const UsdCameraInfo &usdCameraInfo)
{
    pxr::GfVec3f bbMin(std::numeric_limits<float>::max());
    pxr::GfVec3f bbMax(std::numeric_limits<float>::lowest());
    for (const pxr::GfVec3f &p : selected.points)
    {
        bbMin[0] = std::min(bbMin[0], p[0]);
        bbMin[1] = std::min(bbMin[1], p[1]);
        bbMin[2] = std::min(bbMin[2], p[2]);
        bbMax[0] = std::max(bbMax[0], p[0]);
        bbMax[1] = std::max(bbMax[1], p[1]);
        bbMax[2] = std::max(bbMax[2], p[2]);
    }
    const pxr::GfVec3f center = (bbMin + bbMax) * 0.5f;
    const float diag = std::max((bbMax - bbMin).GetLength(), 1.0f);
    const pxr::GfVec3f setupDirection = Normalize(pxr::GfVec3f(0.0f, 0.2f, 2.8f));
    const float fallbackFarDistance = 2.8f * diag;
    const float farDistance =
        usdCameraInfo.found ? std::max(usdCameraInfo.distanceToMeshCenter, 0.001f) : fallbackFarDistance;
    const float nearDistance = farDistance * 0.5f;

    CameraPreset farCamera = {};
    farCamera.name = "far";
    farCamera.origin = center + setupDirection * farDistance;
    farCamera.target = center;
    farCamera.up = pxr::GfVec3f(0.0f, 1.0f, 0.0f);
    farCamera.fovYDegrees = 50.0f;
    farCamera.width = 1920;
    farCamera.height = 1080;
    farCamera.distanceToTarget = farDistance;

    CameraPreset nearCamera = {};
    nearCamera.name = "near";
    nearCamera.origin = center + setupDirection * nearDistance;
    nearCamera.target = center;
    nearCamera.up = pxr::GfVec3f(0.0f, 1.0f, 0.0f);
    nearCamera.fovYDegrees = 50.0f;
    nearCamera.width = 1920;
    nearCamera.height = 1080;
    nearCamera.distanceToTarget = nearDistance;

    return {farCamera, nearCamera};
}

static void BuildPtexMaps(const SelectedSubdivMesh &selected,
                          const Far::PtexIndices &ptexIndices,
                          int numPtexFaces,
                          std::vector<int> &coarseFaceForPtexOut,
                          std::vector<int> &quadrantForPtexOut)
{
    coarseFaceForPtexOut.assign(size_t(numPtexFaces), -1);
    quadrantForPtexOut.assign(size_t(numPtexFaces), -1);

    const int numCoarseFaces = int(selected.faceVertexCounts.size());
    for (int face = 0; face < numCoarseFaces; face++)
    {
        const int basePtex = ptexIndices.GetFaceId(face);
        if (basePtex < 0)
        {
            continue;
        }
        const int count = selected.faceVertexCounts[size_t(face)] == 4
                              ? 1
                              : selected.faceVertexCounts[size_t(face)];
        for (int quadrant = 0; quadrant < count; quadrant++)
        {
            const int ptexFace = basePtex + quadrant;
            if (ptexFace >= 0 && ptexFace < numPtexFaces)
            {
                coarseFaceForPtexOut[size_t(ptexFace)] = face;
                quadrantForPtexOut[size_t(ptexFace)] = quadrant;
            }
        }
    }
}

static std::vector<EdgeAdjacency> BuildAdjacency(const Far::TopologyRefiner &refiner,
                                                 const Far::PtexIndices &ptexIndices,
                                                 const std::vector<int> &coarseFaceForPtex,
                                                 const std::vector<int> &quadrantForPtex)
{
    const int numPtexFaces = int(coarseFaceForPtex.size());
    std::vector<EdgeAdjacency> adjacency(size_t(numPtexFaces), EdgeAdjacency{});
    for (int ptexFace = 0; ptexFace < numPtexFaces; ptexFace++)
    {
        const int coarseFace = coarseFaceForPtex[size_t(ptexFace)];
        const int quadrant = quadrantForPtex[size_t(ptexFace)];
        if (coarseFace < 0 || quadrant < 0)
        {
            continue;
        }

        int adjFaces[4] = {-1, -1, -1, -1};
        int adjEdges[4] = {-1, -1, -1, -1};
        ptexIndices.GetAdjacency(refiner, coarseFace, quadrant, adjFaces, adjEdges);
        for (int edge = 0; edge < 4; edge++)
        {
            adjacency[size_t(ptexFace)].adjacentFace[size_t(edge)] = adjFaces[edge];
            adjacency[size_t(ptexFace)].adjacentEdge[size_t(edge)] = adjEdges[edge];
        }
    }
    return adjacency;
}

static EdgeKey MakeEdgeKey(int face, int edge, int adjacentFace, int adjacentEdge)
{
    if (adjacentFace < 0 || adjacentEdge < 0)
    {
        return {face, edge, -1, -1};
    }

    if (adjacentFace < face || (adjacentFace == face && adjacentEdge < edge))
    {
        return {adjacentFace, adjacentEdge, face, edge};
    }
    return {face, edge, adjacentFace, adjacentEdge};
}

static bool UsesCanonicalOrientation(int face, int edge, const EdgeKey &key)
{
    return face == key.faceA && edge == key.edgeA;
}

static int GetOrCreateEdgeVertex(const EdgeKey &key,
                                 int canonicalSegment,
                                 int segments,
                                 std::unordered_map<EdgeVertexKey, int, EdgeVertexKeyHash> &edgeVertexCache,
                                 TessMesh &mesh,
                                 const Far::PatchMap &patchMap,
                                 const Far::PatchTable &patchTable,
                                 const std::vector<OsdData<pxr::GfVec3f>> &positions)
{
    const EdgeVertexKey vertexKey = {key, canonicalSegment};
    auto it = edgeVertexCache.find(vertexKey);
    if (it != edgeVertexCache.end())
    {
        return it->second;
    }

    const float s = float(canonicalSegment) / float(std::max(1, segments));
    const pxr::GfVec2f uv = EdgeUV(key.edgeA, s);
    const pxr::GfVec3f p = EvaluatePosition(patchMap, patchTable, positions, key.faceA, uv[0], uv[1]);
    const int index = int(mesh.positions.size());
    mesh.positions.push_back(p);
    edgeVertexCache.emplace(vertexKey, index);
    return index;
}

static TessMesh TessellateAdaptiveNoSplit(const Far::PatchMap &patchMap,
                                          const Far::PatchTable &patchTable,
                                          const std::vector<OsdData<pxr::GfVec3f>> &positions,
                                          const std::vector<EdgeAdjacency> &adjacency,
                                          int numPtexFaces,
                                          const CameraPreset &camera,
                                          const AdaptiveSettings &settings)
{
    std::unordered_map<EdgeKey, int, EdgeKeyHash> edgeRateCache;
    std::vector<std::array<int, 4>> patchEdgeRates(
        size_t(numPtexFaces), std::array<int, 4>{{1, 1, 1, 1}});

    for (int face = 0; face < numPtexFaces; face++)
    {
        for (int edge = 0; edge < 4; edge++)
        {
            const int adjFace = adjacency[size_t(face)].adjacentFace[size_t(edge)];
            const int adjEdge = adjacency[size_t(face)].adjacentEdge[size_t(edge)];
            const EdgeKey key = MakeEdgeKey(face, edge, adjFace, adjEdge);
            auto it = edgeRateCache.find(key);
            if (it == edgeRateCache.end())
            {
                const int rate = ComputeEdgeRateTmax(
                    patchMap, patchTable, positions, face, edge, camera, settings);
                it = edgeRateCache.emplace(key, rate).first;
            }
            patchEdgeRates[size_t(face)][size_t(edge)] = it->second;
        }
    }

    TessMesh mesh = {};
    std::unordered_map<EdgeVertexKey, int, EdgeVertexKeyHash> edgeVertexCache;

    for (int face = 0; face < numPtexFaces; face++)
    {
        const std::array<int, 4> &rates = patchEdgeRates[size_t(face)];
        const int uSegments = std::max(rates[0], rates[2]);
        const int vSegments = std::max(rates[1], rates[3]);
        const int width = uSegments + 1;
        const int height = vSegments + 1;

        std::vector<int> gridIndices(size_t(width * height), -1);
        auto GridOffset = [&](int u, int v) -> int { return v * width + u; };

        for (int v = 0; v < height; v++)
        {
            for (int u = 0; u < width; u++)
            {
                const bool onBottom = (v == 0);
                const bool onRight = (u == width - 1);
                const bool onTop = (v == height - 1);
                const bool onLeft = (u == 0);

                if (onBottom || onRight || onTop || onLeft)
                {
                    int edge = 0;
                    float edgeS = 0.0f;
                    if (onBottom)
                    {
                        edge = 0;
                        edgeS = float(u) / float(std::max(1, uSegments));
                    }
                    else if (onRight)
                    {
                        edge = 1;
                        edgeS = float(v) / float(std::max(1, vSegments));
                    }
                    else if (onTop)
                    {
                        edge = 2;
                        edgeS = 1.0f - float(u) / float(std::max(1, uSegments));
                    }
                    else
                    {
                        edge = 3;
                        edgeS = 1.0f - float(v) / float(std::max(1, vSegments));
                    }

                    const int adjFace = adjacency[size_t(face)].adjacentFace[size_t(edge)];
                    const int adjEdge = adjacency[size_t(face)].adjacentEdge[size_t(edge)];
                    const EdgeKey key = MakeEdgeKey(face, edge, adjFace, adjEdge);
                    const int sharedRate = patchEdgeRates[size_t(face)][size_t(edge)];
                    const int localSegment =
                        std::clamp(int(std::lround(edgeS * float(sharedRate))), 0, sharedRate);
                    const int canonicalSegment = UsesCanonicalOrientation(face, edge, key)
                                                     ? localSegment
                                                     : (sharedRate - localSegment);

                    const int index = GetOrCreateEdgeVertex(key,
                                                            canonicalSegment,
                                                            sharedRate,
                                                            edgeVertexCache,
                                                            mesh,
                                                            patchMap,
                                                            patchTable,
                                                            positions);
                    gridIndices[size_t(GridOffset(u, v))] = index;
                }
                else
                {
                    const float fu = float(u) / float(std::max(1, uSegments));
                    const float fv = float(v) / float(std::max(1, vSegments));
                    const pxr::GfVec3f p =
                        EvaluatePosition(patchMap, patchTable, positions, face, fu, fv);
                    const int index = int(mesh.positions.size());
                    mesh.positions.push_back(p);
                    gridIndices[size_t(GridOffset(u, v))] = index;
                }
            }
        }

        for (int v = 0; v < vSegments; v++)
        {
            for (int u = 0; u < uSegments; u++)
            {
                const int i0 = gridIndices[size_t(GridOffset(u, v))];
                const int i1 = gridIndices[size_t(GridOffset(u + 1, v))];
                const int i2 = gridIndices[size_t(GridOffset(u + 1, v + 1))];
                const int i3 = gridIndices[size_t(GridOffset(u, v + 1))];
                mesh.indices.push_back(i0);
                mesh.indices.push_back(i1);
                mesh.indices.push_back(i2);
                mesh.indices.push_back(i0);
                mesh.indices.push_back(i2);
                mesh.indices.push_back(i3);
            }
        }
    }

    return mesh;
}

static bool WriteAdaptiveObj(const TessMesh &mesh,
                             const fs::path &path,
                             const SelectedSubdivMesh &source,
                             const CameraPreset &camera,
                             const AdaptiveSettings &settings,
                             const UsdCameraInfo &usdCameraInfo)
{
    std::ofstream out(path, std::ios::out | std::ios::binary);
    if (!out.is_open())
    {
        return false;
    }

    out << "# source_prim " << source.path.GetString() << "\n";
    out << "# scheme " << source.subdivisionScheme << "\n";
    out << "# control_cage_faces " << source.faceVertexCounts.size() << "\n";
    out << "# mode adaptive_nosplit\n";
    out << "# camera " << camera.name << "\n";
    out << "# camera_distance_to_target " << camera.distanceToTarget << "\n";
    if (usdCameraInfo.found)
    {
        out << "# usd_camera_path " << usdCameraInfo.path.GetString() << "\n";
        out << "# usd_camera_distance_to_mesh_center " << usdCameraInfo.distanceToMeshCenter << "\n";
    }
    out << "# N " << settings.numEdgeSamples << "\n";
    out << "# R " << settings.pixelSpacing << "\n";
    out << "# rate_clamp_min " << settings.minRate << "\n";
    out << "# rate_clamp_max " << settings.maxRate << "\n";
    out << "# vertices " << mesh.positions.size() << "\n";
    out << "# triangles " << (mesh.indices.size() / 3) << "\n";

    for (const pxr::GfVec3f &p : mesh.positions)
    {
        out << "v " << p[0] << " " << p[1] << " " << p[2] << "\n";
    }
    for (size_t i = 0; i + 2 < mesh.indices.size(); i += 3)
    {
        out << "f " << (mesh.indices[i + 0] + 1) << " " << (mesh.indices[i + 1] + 1) << " "
            << (mesh.indices[i + 2] + 1) << "\n";
    }
    return out.good();
}

static bool WriteControlCageObj(const SelectedSubdivMesh &source, const fs::path &path)
{
    std::ofstream out(path, std::ios::out | std::ios::binary);
    if (!out.is_open())
    {
        return false;
    }

    out << "# source_prim " << source.path.GetString() << "\n";
    out << "# scheme " << source.subdivisionScheme << "\n";
    out << "# control_cage_faces " << source.faceVertexCounts.size() << "\n";
    out << "# control_cage_vertices " << source.points.size() << "\n";

    for (const pxr::GfVec3f &p : source.points)
    {
        out << "v " << p[0] << " " << p[1] << " " << p[2] << "\n";
    }

    int indexOffset = 0;
    for (int faceVertexCount : source.faceVertexCounts)
    {
        if (faceVertexCount < 3)
        {
            indexOffset += faceVertexCount;
            continue;
        }

        const int i0 = source.faceVertexIndices[indexOffset];
        for (int i = 1; i < faceVertexCount - 1; i++)
        {
            const int i1 = source.faceVertexIndices[indexOffset + i];
            const int i2 = source.faceVertexIndices[indexOffset + i + 1];
            out << "f " << (i0 + 1) << " " << (i1 + 1) << " " << (i2 + 1) << "\n";
        }
        indexOffset += faceVertexCount;
    }

    return out.good();
}

static bool ProcessSelectedMesh(const pxr::UsdStageRefPtr &stage,
                                const SelectedSubdivMesh &selected,
                                const fs::path &outputDir,
                                const char *controlCageName,
                                const char *farOutputName,
                                const char *nearOutputName)
{
    Far::TopologyDescriptor desc = {};
    desc.numVertices = int(selected.points.size());
    desc.numFaces = int(selected.faceVertexCounts.size());
    desc.numVertsPerFace = selected.faceVertexCounts.cdata();
    desc.vertIndicesPerFace = selected.faceVertexIndices.cdata();
    std::vector<int> creaseVertexPairs;
    std::vector<float> creaseSharpness;
    BuildCreasePairs(selected, creaseVertexPairs, creaseSharpness);
    desc.numCreases = int(creaseSharpness.size());
    desc.creaseVertexIndexPairs = creaseVertexPairs.empty() ? nullptr : creaseVertexPairs.data();
    desc.creaseWeights = creaseSharpness.empty() ? nullptr : creaseSharpness.data();
    desc.numCorners = int(std::min(selected.cornerIndices.size(), selected.cornerSharpnesses.size()));
    desc.cornerVertexIndices = selected.cornerIndices.cdata();
    desc.cornerWeights = selected.cornerSharpnesses.cdata();
    desc.numHoles = int(selected.holeIndices.size());
    desc.holeIndices = selected.holeIndices.cdata();

    const Sdc::SchemeType scheme = Sdc::SCHEME_CATMARK;
    Sdc::Options options;
    Far::TopologyRefiner *refiner = Far::TopologyRefinerFactory<Far::TopologyDescriptor>::Create(
        desc, Far::TopologyRefinerFactory<Far::TopologyDescriptor>::Options(scheme, options));
    if (!refiner)
    {
        fprintf(stderr, "Failed to create topology refiner.\n");
        return false;
    }

    Far::TopologyRefiner::AdaptiveOptions adaptiveOptions(1);
    refiner->RefineAdaptive(adaptiveOptions);

    Far::PatchTableFactory::Options patchOptions(1);
    patchOptions.endCapType = Far::PatchTableFactory::Options::ENDCAP_GREGORY_BASIS;
    const Far::PatchTable *patchTable = Far::PatchTableFactory::Create(*refiner, patchOptions);
    if (!patchTable)
    {
        fprintf(stderr, "Failed to create patch table.\n");
        delete refiner;
        return false;
    }

    Far::PatchMap patchMap(*patchTable);
    Far::PtexIndices ptexIndices(*refiner);
    const int numPtexFaces = ptexIndices.GetNumFaces();

    std::vector<int> coarseFaceForPtex;
    std::vector<int> quadrantForPtex;
    BuildPtexMaps(selected, ptexIndices, numPtexFaces, coarseFaceForPtex, quadrantForPtex);
    const std::vector<EdgeAdjacency> adjacency =
        BuildAdjacency(*refiner, ptexIndices, coarseFaceForPtex, quadrantForPtex);

    std::vector<pxr::GfVec3f> cagePositions;
    cagePositions.reserve(selected.points.size());
    for (const pxr::GfVec3f &p : selected.points)
    {
        cagePositions.push_back(p);
    }
    const std::vector<OsdData<pxr::GfVec3f>> refinedPositions =
        InterpolateVertex(cagePositions, refiner, patchTable);

    const fs::path controlCagePath = outputDir / controlCageName;
    if (!WriteControlCageObj(selected, controlCagePath))
    {
        fprintf(stderr, "Failed to write control cage OBJ: %s\n", controlCagePath.string().c_str());
        delete patchTable;
        delete refiner;
        return false;
    }
    printf("Wrote %s (verts=%zu faces=%zu)\n",
           controlCagePath.string().c_str(),
           selected.points.size(),
           selected.faceVertexCounts.size());
    printf("Control-cage features: creaseSegments=%zu corners=%d holes=%d\n",
           creaseSharpness.size(),
           desc.numCorners,
           desc.numHoles);

    UsdCameraInfo usdCameraInfo = {};
    if (GetClosestUsdCameraInfo(stage, selected.points, usdCameraInfo))
    {
        printf("Using USD camera: %s distance=%f\n",
               usdCameraInfo.path.GetString().c_str(),
               usdCameraInfo.distanceToMeshCenter);
    }
    else
    {
        printf("No USD camera found; using fallback mesh-scaled camera distance.\n");
    }

    const std::array<CameraPreset, 2> cameras = BuildCameraPresets(selected, usdCameraInfo);
    const AdaptiveSettings settings = {};
    const char *outputNames[2] = {farOutputName, nearOutputName};
    for (int i = 0; i < 2; i++)
    {
        const CameraPreset &camera = cameras[size_t(i)];
        const TessMesh tessMesh = TessellateAdaptiveNoSplit(
            patchMap, *patchTable, refinedPositions, adjacency, numPtexFaces, camera, settings);
        const fs::path outPath = outputDir / outputNames[i];
        if (!WriteAdaptiveObj(tessMesh, outPath, selected, camera, settings, usdCameraInfo))
        {
            fprintf(stderr, "Failed to write OBJ: %s\n", outPath.string().c_str());
            delete patchTable;
            delete refiner;
            return false;
        }
        printf("Wrote %s (camera=%s verts=%zu tris=%zu)\n",
               outPath.string().c_str(),
               camera.name,
               tessMesh.positions.size(),
               tessMesh.indices.size() / 3);
    }

    delete patchTable;
    delete refiner;
    return true;
}

int main()
{
    const std::string usdPath = "C:/Users/maven/Downloads/ALab-2.2.0/ALab/entry.usda";
    pxr::UsdStageRefPtr stage = pxr::UsdStage::Open(usdPath);
    if (!stage)
    {
        fprintf(stderr, "Failed to open USD stage: %s\n", usdPath.c_str());
        return 1;
    }

    const fs::path outputDir = fs::path("tests") / "bvh" / "out";
    fs::create_directories(outputDir);

    SelectedSubdivMesh selected = {};
    if (!SelectLargestCatmullClarkMesh(stage, selected))
    {
        fprintf(stderr, "No Catmull-Clark mesh found in stage.\n");
        return 1;
    }
    if (!ProcessSelectedMesh(stage,
                             selected,
                             outputDir,
                             "subdiv_control_cage.obj",
                             "subdiv_limit_adaptive_nosplit_far.obj",
                             "subdiv_limit_adaptive_nosplit_near.obj"))
    {
        return 1;
    }

    SelectedSubdivMesh selectedWithCreases = {};
    if (SelectLargestCatmullClarkMeshWithCreases(stage, selectedWithCreases))
    {
        if (!ProcessSelectedMesh(stage,
                                 selectedWithCreases,
                                 outputDir,
                                 "subdiv_control_cage_creased.obj",
                                 "subdiv_limit_adaptive_nosplit_far_creased.obj",
                                 "subdiv_limit_adaptive_nosplit_near_creased.obj"))
        {
            return 1;
        }
    }
    else
    {
        printf("No Catmull-Clark mesh with creases found in stage.\n");
    }

    return 0;
}
