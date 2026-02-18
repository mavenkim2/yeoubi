#include "bvh/usd_subdiv_select.h"
#include "bvh/usd_camera_utils.h"
#include "io/usd_mesh_io.h"
#include "io/usd_subdiv_json_io.h"

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
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <numeric>
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

static const SelectedSubdivMesh *gFallbackSelectedMesh = nullptr;
static const std::vector<int> *gFallbackCoarseFaceForPtex = nullptr;
static const std::vector<int> *gFallbackFaceStartOffsets = nullptr;

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

enum class EvaluationMode
{
    Adaptive,
    Uniform
};

struct CliOptions
{
    fs::path inputJsonPath;
    fs::path outputDir = fs::path("tests") / "bvh" / "out";
    std::string farOutputName = "subdiv_limit_adaptive_nosplit_far.obj";
    std::string nearOutputName = "subdiv_limit_adaptive_nosplit_near.obj";
    EvaluationMode evaluationMode = EvaluationMode::Adaptive;
    int uniformRate = 8;
};

static void PrintUsage(const char *exeName)
{
    printf("Usage: %s --input-json <path> [--out-dir <dir>] [--far-out <obj>] [--near-out <obj>] "
           "[--evaluation adaptive|uniform] [--uniform-rate N]\n",
           exeName);
}

static CliOptions ParseCli(int argc, char **argv)
{
    CliOptions options = {};
    for (int i = 1; i < argc; i++)
    {
        const std::string arg = argv[i];
        if (arg == "--input-json")
        {
            if (i + 1 >= argc)
            {
                PrintUsage(argv[0]);
                std::abort();
            }
            options.inputJsonPath = fs::path(argv[++i]);
            continue;
        }
        if (arg == "--out-dir")
        {
            if (i + 1 >= argc)
            {
                PrintUsage(argv[0]);
                std::abort();
            }
            options.outputDir = fs::path(argv[++i]);
            continue;
        }
        if (arg == "--far-out")
        {
            if (i + 1 >= argc)
            {
                PrintUsage(argv[0]);
                std::abort();
            }
            options.farOutputName = argv[++i];
            continue;
        }
        if (arg == "--near-out")
        {
            if (i + 1 >= argc)
            {
                PrintUsage(argv[0]);
                std::abort();
            }
            options.nearOutputName = argv[++i];
            continue;
        }
        if (arg == "--evaluation")
        {
            if (i + 1 >= argc)
            {
                PrintUsage(argv[0]);
                std::abort();
            }
            const std::string mode = argv[++i];
            if (mode == "adaptive")
            {
                options.evaluationMode = EvaluationMode::Adaptive;
            }
            else if (mode == "uniform")
            {
                options.evaluationMode = EvaluationMode::Uniform;
            }
            else
            {
                PrintUsage(argv[0]);
                std::abort();
            }
            continue;
        }
        if (arg == "--uniform-rate")
        {
            if (i + 1 >= argc)
            {
                PrintUsage(argv[0]);
                std::abort();
            }
            options.uniformRate = std::max(1, std::stoi(argv[++i]));
            continue;
        }
        if (arg == "--help" || arg == "-h")
        {
            PrintUsage(argv[0]);
            std::exit(0);
        }

        PrintUsage(argv[0]);
        std::abort();
    }

    if (options.inputJsonPath.empty())
    {
        PrintUsage(argv[0]);
        std::abort();
    }
    return options;
}


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

static const Far::PatchTable::PatchHandle *FindPatchHandleRobust(const Far::PatchMap &patchMap,
                                                                  int ptexFace,
                                                                  float u,
                                                                  float v,
                                                                  float &evalUOut,
                                                                  float &evalVOut)
{
    const float eps = 1e-6f;
    const float clampedU = std::clamp(u, eps, 1.0f - eps);
    const float clampedV = std::clamp(v, eps, 1.0f - eps);
    const std::array<pxr::GfVec2f, 12> candidates = {
        pxr::GfVec2f(u, v),
        pxr::GfVec2f(clampedU, clampedV),
        pxr::GfVec2f(std::clamp(clampedU + eps, eps, 1.0f - eps), clampedV),
        pxr::GfVec2f(std::clamp(clampedU - eps, eps, 1.0f - eps), clampedV),
        pxr::GfVec2f(clampedU, std::clamp(clampedV + eps, eps, 1.0f - eps)),
        pxr::GfVec2f(clampedU, std::clamp(clampedV - eps, eps, 1.0f - eps)),
        pxr::GfVec2f(std::clamp(clampedU + eps, eps, 1.0f - eps),
                     std::clamp(clampedV + eps, eps, 1.0f - eps)),
        pxr::GfVec2f(std::clamp(clampedU + eps, eps, 1.0f - eps),
                     std::clamp(clampedV - eps, eps, 1.0f - eps)),
        pxr::GfVec2f(std::clamp(clampedU - eps, eps, 1.0f - eps),
                     std::clamp(clampedV + eps, eps, 1.0f - eps)),
        pxr::GfVec2f(std::clamp(clampedU - eps, eps, 1.0f - eps),
                     std::clamp(clampedV - eps, eps, 1.0f - eps)),
        pxr::GfVec2f(0.5f, 0.5f),
        pxr::GfVec2f(0.25f, 0.25f),
    };

    for (const pxr::GfVec2f &candidate : candidates)
    {
        const Far::PatchTable::PatchHandle *handle =
            patchMap.FindPatch(ptexFace, candidate[0], candidate[1]);
        if (handle)
        {
            evalUOut = candidate[0];
            evalVOut = candidate[1];
            return handle;
        }
    }
    evalUOut = clampedU;
    evalVOut = clampedV;
    return nullptr;
}

static pxr::GfVec3f EvaluateCoarseFallbackPosition(int ptexFace, float u, float v)
{
    if (!gFallbackSelectedMesh || !gFallbackCoarseFaceForPtex || !gFallbackFaceStartOffsets ||
        ptexFace < 0 || ptexFace >= int(gFallbackCoarseFaceForPtex->size()))
    {
        return pxr::GfVec3f(0.0f);
    }

    const int coarseFace = (*gFallbackCoarseFaceForPtex)[size_t(ptexFace)];
    if (coarseFace < 0 || coarseFace >= int(gFallbackSelectedMesh->faceVertexCounts.size()))
    {
        return pxr::GfVec3f(0.0f);
    }

    const int faceVertexCount = gFallbackSelectedMesh->faceVertexCounts[size_t(coarseFace)];
    if (faceVertexCount != 4)
    {
        return pxr::GfVec3f(0.0f);
    }

    const int faceStart = (*gFallbackFaceStartOffsets)[size_t(coarseFace)];
    const int i0 = gFallbackSelectedMesh->faceVertexIndices[size_t(faceStart + 0)];
    const int i1 = gFallbackSelectedMesh->faceVertexIndices[size_t(faceStart + 1)];
    const int i2 = gFallbackSelectedMesh->faceVertexIndices[size_t(faceStart + 2)];
    const int i3 = gFallbackSelectedMesh->faceVertexIndices[size_t(faceStart + 3)];
    const pxr::GfVec3f p0 = gFallbackSelectedMesh->points[size_t(i0)];
    const pxr::GfVec3f p1 = gFallbackSelectedMesh->points[size_t(i1)];
    const pxr::GfVec3f p2 = gFallbackSelectedMesh->points[size_t(i2)];
    const pxr::GfVec3f p3 = gFallbackSelectedMesh->points[size_t(i3)];
    const float uu = std::clamp(u, 0.0f, 1.0f);
    const float vv = std::clamp(v, 0.0f, 1.0f);
    const pxr::GfVec3f bottom = p0 * (1.0f - uu) + p1 * uu;
    const pxr::GfVec3f top = p3 * (1.0f - uu) + p2 * uu;
    return bottom * (1.0f - vv) + top * vv;
}

static pxr::GfVec3f EvaluatePosition(const Far::PatchMap &patchMap,
                                     const Far::PatchTable &patchTable,
                                     const std::vector<OsdData<pxr::GfVec3f>> &positions,
                                     int ptexFace,
                                     float u,
                                     float v)
{
    const Far::PatchTable::PatchHandle *rawHandle = patchMap.FindPatch(ptexFace, u, v);
    float evalU = 0.0f;
    float evalV = 0.0f;
    const Far::PatchTable::PatchHandle *handle =
        FindPatchHandleRobust(patchMap, ptexFace, u, v, evalU, evalV);
    (void)rawHandle;
    if (!handle)
    {
        const pxr::GfVec3f coarseFallback = EvaluateCoarseFallbackPosition(ptexFace, u, v);
        if (!(fabsf(coarseFallback[0]) < 1e-20f && fabsf(coarseFallback[1]) < 1e-20f &&
              fabsf(coarseFallback[2]) < 1e-20f))
        {
            return coarseFallback;
        }
        if (!positions.empty())
        {
            return positions[0].value;
        }
        return pxr::GfVec3f(0.0f);
    }

    Far::ConstIndexArray cvIndices = patchTable.GetPatchVertices(*handle);
    std::vector<float> pWeights(size_t(cvIndices.size()));
    patchTable.EvaluateBasis(*handle, evalU, evalV, pWeights.data());

    pxr::GfVec3f p(0.0f);
    for (int cv = 0; cv < cvIndices.size(); cv++)
    {
        p += positions[size_t(cvIndices[cv])].value * pWeights[size_t(cv)];
    }
    return p;
}

static bool FaceHasPatchCoverage(const Far::PatchMap &patchMap, int ptexFace)
{
    const std::array<pxr::GfVec2f, 5> probes = {
        pxr::GfVec2f(0.5f, 0.5f),
        pxr::GfVec2f(0.25f, 0.25f),
        pxr::GfVec2f(0.75f, 0.25f),
        pxr::GfVec2f(0.25f, 0.75f),
        pxr::GfVec2f(0.75f, 0.75f)};
    for (const pxr::GfVec2f &probe : probes)
    {
        float evalU = 0.0f;
        float evalV = 0.0f;
        if (FindPatchHandleRobust(patchMap, ptexFace, probe[0], probe[1], evalU, evalV))
        {
            return true;
        }
    }
    return false;
}

static void PrintNonEvaluableFaceInfo(const Far::PatchMap &patchMap,
                                      const std::vector<int> &coarseFaceForPtex,
                                      const SelectedSubdivMesh &selected,
                                      const char *label)
{
    std::unordered_map<int, int> nonEvaluableCountByFaceVertexCount;
    int totalNonEvaluable = 0;
    for (size_t ptex = 0; ptex < coarseFaceForPtex.size(); ptex++)
    {
        if (FaceHasPatchCoverage(patchMap, int(ptex)))
        {
            continue;
        }

        totalNonEvaluable++;
        const int coarseFace = coarseFaceForPtex[ptex];
        int faceVertexCount = -1;
        if (coarseFace >= 0 && coarseFace < int(selected.faceVertexCounts.size()))
        {
            faceVertexCount = selected.faceVertexCounts[size_t(coarseFace)];
        }
        nonEvaluableCountByFaceVertexCount[faceVertexCount]++;
        if (totalNonEvaluable <= 128)
        {
            printf("non_evaluable [%s] ptex=%zu coarse_face=%d face_vertex_count=%d\n",
                   label,
                   ptex,
                   coarseFace,
                   faceVertexCount);
        }
    }

    printf("non_evaluable_summary [%s] total=%d\n", label, totalNonEvaluable);
    for (const auto &entry : nonEvaluableCountByFaceVertexCount)
    {
        printf("  non_evaluable_face_vertex_count [%s] count=%d faces=%d\n",
               label,
               entry.first,
               entry.second);
    }
}

static void PrintRawPatchLookupProbe(const Far::PatchMap &patchMap,
                                     int numPtexFaces,
                                     const char *label)
{
    int missCount = 0;
    int printed = 0;
    for (int face = 0; face < numPtexFaces; face++)
    {
        for (int vStep = 0; vStep <= 8; vStep++)
        {
            for (int uStep = 0; uStep <= 8; uStep++)
            {
                const float u = float(uStep) / 8.0f;
                const float v = float(vStep) / 8.0f;
                const Far::PatchTable::PatchHandle *handle = patchMap.FindPatch(face, u, v);
                if (!handle)
                {
                    missCount++;
                    if (printed < 24)
                    {
                        printf("raw_lookup_miss [%s] ptex=%d uv=(%.3f,%.3f)\n",
                               label,
                               face,
                               u,
                               v);
                        printed++;
                    }
                }
            }
        }
    }
    printf("raw_lookup_probe [%s]: misses=%d sampled=%d\n",
           label,
           missCount,
           numPtexFaces * 81);
}

static bool FaceHasAnyRawPatchCoverageDense(const Far::PatchMap &patchMap, int ptexFace)
{
    for (int vStep = 0; vStep <= 8; vStep++)
    {
        for (int uStep = 0; uStep <= 8; uStep++)
        {
            const float u = float(uStep) / 8.0f;
            const float v = float(vStep) / 8.0f;
            if (patchMap.FindPatch(ptexFace, u, v))
            {
                return true;
            }
        }
    }
    return false;
}

static void PrintRawNonEvaluableFaceInfo(const Far::PatchMap &patchMap,
                                         const std::vector<int> &coarseFaceForPtex,
                                         const SelectedSubdivMesh &selected,
                                         const char *label)
{
    if (selected.faceVertexCounts.size() > 417)
    {
        const int faceVertexCount417 = selected.faceVertexCounts[417];
        printf("coarse_face_probe [%s] coarse_face=417 face_vertex_count=%d is_triangle=%d\n",
               label,
               faceVertexCount417,
               faceVertexCount417 == 3 ? 1 : 0);
    }

    std::unordered_map<int, int> nonEvaluableCountByFaceVertexCount;
    int totalRawNonEvaluable = 0;
    for (size_t ptex = 0; ptex < coarseFaceForPtex.size(); ptex++)
    {
        if (FaceHasAnyRawPatchCoverageDense(patchMap, int(ptex)))
        {
            continue;
        }

        totalRawNonEvaluable++;
        const int coarseFace = coarseFaceForPtex[ptex];
        int faceVertexCount = -1;
        if (coarseFace >= 0 && coarseFace < int(selected.faceVertexCounts.size()))
        {
            faceVertexCount = selected.faceVertexCounts[size_t(coarseFace)];
        }
        nonEvaluableCountByFaceVertexCount[faceVertexCount]++;
        printf("raw_non_evaluable [%s] ptex=%zu coarse_face=%d face_vertex_count=%d\n",
               label,
               ptex,
               coarseFace,
               faceVertexCount);
    }

    printf("raw_non_evaluable_summary [%s] total=%d\n", label, totalRawNonEvaluable);
    for (const auto &entry : nonEvaluableCountByFaceVertexCount)
    {
        printf("  raw_non_evaluable_face_vertex_count [%s] count=%d faces=%d\n",
               label,
               entry.first,
               entry.second);
    }

    if (420 < int(coarseFaceForPtex.size()))
    {
        const int coarseFace420 = coarseFaceForPtex[420];
        int faceVertexCount420 = -1;
        if (coarseFace420 >= 0 && coarseFace420 < int(selected.faceVertexCounts.size()))
        {
            faceVertexCount420 = selected.faceVertexCounts[size_t(coarseFace420)];
        }
        printf("raw_non_evaluable_probe ptex=420 [%s] coarse_face=%d face_vertex_count=%d has_raw_coverage=%d\n",
               label,
               coarseFace420,
               faceVertexCount420,
               FaceHasAnyRawPatchCoverageDense(patchMap, 420) ? 1 : 0);
    }
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
                                 const std::vector<bool> &faceIsEvaluable,
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

    int evalFace = key.faceA;
    int evalEdge = key.edgeA;
    int evalSegment = canonicalSegment;
    if (key.faceB >= 0 && key.edgeB >= 0)
    {
        const bool faceAEvaluable =
            key.faceA >= 0 && key.faceA < int(faceIsEvaluable.size()) && faceIsEvaluable[size_t(key.faceA)];
        const bool faceBEvaluable =
            key.faceB >= 0 && key.faceB < int(faceIsEvaluable.size()) && faceIsEvaluable[size_t(key.faceB)];
        if (!faceAEvaluable && faceBEvaluable)
        {
            evalFace = key.faceB;
            evalEdge = key.edgeB;
            evalSegment = segments - canonicalSegment;
        }
    }

    const float s = float(evalSegment) / float(std::max(1, segments));
    const pxr::GfVec2f uv = EdgeUV(evalEdge, s);
    const pxr::GfVec3f p = EvaluatePosition(patchMap, patchTable, positions, evalFace, uv[0], uv[1]);
    const int index = int(mesh.positions.size());
    mesh.positions.push_back(p);
    edgeVertexCache.emplace(vertexKey, index);
    return index;
}

static float TriangleAreaSquared(const pxr::GfVec3f &a, const pxr::GfVec3f &b, const pxr::GfVec3f &c)
{
    const pxr::GfVec3f ab = b - a;
    const pxr::GfVec3f ac = c - a;
    const pxr::GfVec3f cross = pxr::GfCross(ab, ac);
    return pxr::GfDot(cross, cross) * 0.25f;
}

static bool IsTriangleDegenerate(const TessMesh &mesh, int i0, int i1, int i2)
{
    if (i0 == i1 || i1 == i2 || i2 == i0)
    {
        return true;
    }
    if (i0 < 0 || i1 < 0 || i2 < 0 || i0 >= int(mesh.positions.size()) || i1 >= int(mesh.positions.size()) ||
        i2 >= int(mesh.positions.size()))
    {
        return true;
    }

    const float areaSq = TriangleAreaSquared(
        mesh.positions[size_t(i0)], mesh.positions[size_t(i1)], mesh.positions[size_t(i2)]);
    return areaSq < 1e-18f;
}

static void EmitQuadCellTriangles(TessMesh &mesh, int i0, int i1, int i2, int i3)
{
    std::array<int, 4> ring = {i0, i1, i2, i3};
    std::vector<int> unique;
    unique.reserve(4);
    for (int i = 0; i < 4; i++)
    {
        const int index = ring[size_t(i)];
        if (unique.empty() || unique.back() != index)
        {
            unique.push_back(index);
        }
    }
    if (unique.size() >= 2 && unique.front() == unique.back())
    {
        unique.pop_back();
    }

    if (unique.size() < 3)
    {
        return;
    }
    if (unique.size() == 3)
    {
        if (!IsTriangleDegenerate(mesh, unique[0], unique[1], unique[2]))
        {
            mesh.indices.push_back(unique[0]);
            mesh.indices.push_back(unique[1]);
            mesh.indices.push_back(unique[2]);
        }
        return;
    }

    if (!IsTriangleDegenerate(mesh, unique[0], unique[1], unique[2]))
    {
        mesh.indices.push_back(unique[0]);
        mesh.indices.push_back(unique[1]);
        mesh.indices.push_back(unique[2]);
    }
    if (!IsTriangleDegenerate(mesh, unique[0], unique[2], unique[3]))
    {
        mesh.indices.push_back(unique[0]);
        mesh.indices.push_back(unique[2]);
        mesh.indices.push_back(unique[3]);
    }
}

static void EmitTriangleIfValid(TessMesh &mesh, int i0, int i1, int i2)
{
    if (!IsTriangleDegenerate(mesh, i0, i1, i2))
    {
        mesh.indices.push_back(i0);
        mesh.indices.push_back(i1);
        mesh.indices.push_back(i2);
    }
}

static void EmitLowEdgeStitchTriangles(TessMesh &mesh,
                                       const std::vector<int> &gridIndices,
                                       int uSegments,
                                       int vSegments,
                                       int edge,
                                       int lowRate)
{
    auto GridOffset = [&](int u, int v) -> int { return v * (uSegments + 1) + u; };
    auto BoundaryVertex = [&](int t) -> int {
        if (edge == 0)
        {
            return gridIndices[size_t(GridOffset(t, 0))];
        }
        if (edge == 1)
        {
            return gridIndices[size_t(GridOffset(uSegments, t))];
        }
        if (edge == 2)
        {
            return gridIndices[size_t(GridOffset(uSegments - t, vSegments))];
        }
        return gridIndices[size_t(GridOffset(0, vSegments - t))];
    };
    auto InnerVertex = [&](int t) -> int {
        if (edge == 0)
        {
            return gridIndices[size_t(GridOffset(t, 1))];
        }
        if (edge == 1)
        {
            return gridIndices[size_t(GridOffset(uSegments - 1, t))];
        }
        if (edge == 2)
        {
            return gridIndices[size_t(GridOffset(uSegments - t, vSegments - 1))];
        }
        return gridIndices[size_t(GridOffset(1, vSegments - t))];
    };

    const int highRate = (edge == 0 || edge == 2) ? uSegments : vSegments;
    if (highRate <= 0 || lowRate <= 0)
    {
        return;
    }

    for (int coarseSegment = 0; coarseSegment < lowRate; coarseSegment++)
    {
        int h0 = int(std::lround(float(coarseSegment) * float(highRate) / float(lowRate)));
        int h1 = int(std::lround(float(coarseSegment + 1) * float(highRate) / float(lowRate)));
        h0 = std::clamp(h0, 0, highRate);
        h1 = std::clamp(h1, 0, highRate);
        if (h1 <= h0)
        {
            h1 = std::min(highRate, h0 + 1);
        }
        const int tip0 = BoundaryVertex(h0);
        const int tip1 = BoundaryVertex(h1);
        for (int h = h0; h < h1; h++)
        {
            const int i0 = InnerVertex(h);
            const int i1 = InnerVertex(h + 1);
            EmitTriangleIfValid(mesh, tip0, i1, i0);
        }
        if (tip1 != tip0)
        {
            const int edgeInnerEnd = InnerVertex(h1);
            EmitTriangleIfValid(mesh, tip0, tip1, edgeInnerEnd);
        }
    }
}

static TessMesh TessellateAdaptiveNoSplit(const Far::PatchMap &patchMap,
                                          const Far::PatchTable &patchTable,
                                          const std::vector<OsdData<pxr::GfVec3f>> &positions,
                                          const std::vector<EdgeAdjacency> &adjacency,
                                          const std::vector<int> &coarseFaceForPtex,
                                          const std::vector<int> &faceStartOffsets,
                                          const SelectedSubdivMesh &selected,
                                          int numPtexFaces,
                                          const CameraPreset &camera,
                                          const AdaptiveSettings &settings)
{
    std::vector<bool> faceIsEvaluable(size_t(numPtexFaces), false);
    int nonEvaluableFaceCount = 0;
    for (int face = 0; face < numPtexFaces; face++)
    {
        faceIsEvaluable[size_t(face)] = FaceHasPatchCoverage(patchMap, face);
        if (!faceIsEvaluable[size_t(face)])
        {
            nonEvaluableFaceCount++;
        }
    }
    if (nonEvaluableFaceCount > 0)
    {
        printf("Tessellation warning: %d/%d ptex faces are non-evaluable; using coarse fallback for non-evaluable samples.\n",
               nonEvaluableFaceCount,
               numPtexFaces);
    }
    gFallbackSelectedMesh = &selected;
    gFallbackCoarseFaceForPtex = &coarseFaceForPtex;
    gFallbackFaceStartOffsets = &faceStartOffsets;

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
        const bool lowBottom = rates[0] < rates[2];
        const bool lowRight = rates[1] < rates[3];
        const bool lowTop = rates[2] < rates[0];
        const bool lowLeft = rates[3] < rates[1];
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
                                                            faceIsEvaluable,
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
                if (lowBottom && v == 0)
                {
                    continue;
                }
                if (lowTop && v == vSegments - 1)
                {
                    continue;
                }
                if (lowLeft && u == 0)
                {
                    continue;
                }
                if (lowRight && u == uSegments - 1)
                {
                    continue;
                }
                const int i0 = gridIndices[size_t(GridOffset(u, v))];
                const int i1 = gridIndices[size_t(GridOffset(u + 1, v))];
                const int i2 = gridIndices[size_t(GridOffset(u + 1, v + 1))];
                const int i3 = gridIndices[size_t(GridOffset(u, v + 1))];
                EmitQuadCellTriangles(mesh, i0, i1, i2, i3);
            }
        }
        if (lowBottom)
        {
            EmitLowEdgeStitchTriangles(mesh, gridIndices, uSegments, vSegments, 0, rates[0]);
        }
        if (lowRight)
        {
            EmitLowEdgeStitchTriangles(mesh, gridIndices, uSegments, vSegments, 1, rates[1]);
        }
        if (lowTop)
        {
            EmitLowEdgeStitchTriangles(mesh, gridIndices, uSegments, vSegments, 2, rates[2]);
        }
        if (lowLeft)
        {
            EmitLowEdgeStitchTriangles(mesh, gridIndices, uSegments, vSegments, 3, rates[3]);
        }
    }

    gFallbackSelectedMesh = nullptr;
    gFallbackCoarseFaceForPtex = nullptr;
    gFallbackFaceStartOffsets = nullptr;

    return mesh;
}

static TessMesh TessellateUniformNoStitch(const Far::PatchMap &patchMap,
                                          const Far::PatchTable &patchTable,
                                          const std::vector<OsdData<pxr::GfVec3f>> &positions,
                                          const std::vector<int> &coarseFaceForPtex,
                                          const std::vector<int> &faceStartOffsets,
                                          const SelectedSubdivMesh &selected,
                                          int numPtexFaces,
                                          int uniformRate)
{
    gFallbackSelectedMesh = &selected;
    gFallbackCoarseFaceForPtex = &coarseFaceForPtex;
    gFallbackFaceStartOffsets = &faceStartOffsets;

    TessMesh mesh = {};
    const int segments = std::max(1, uniformRate);
    const int width = segments + 1;
    const int height = segments + 1;
    for (int face = 0; face < numPtexFaces; face++)
    {
        std::vector<int> gridIndices(size_t(width * height), -1);
        auto GridOffset = [&](int u, int v) -> int { return v * width + u; };
        for (int v = 0; v < height; v++)
        {
            for (int u = 0; u < width; u++)
            {
                const float fu = float(u) / float(segments);
                const float fv = float(v) / float(segments);
                const pxr::GfVec3f p = EvaluatePosition(patchMap, patchTable, positions, face, fu, fv);
                const int index = int(mesh.positions.size());
                mesh.positions.push_back(p);
                gridIndices[size_t(GridOffset(u, v))] = index;
            }
        }

        for (int v = 0; v < segments; v++)
        {
            for (int u = 0; u < segments; u++)
            {
                const int i0 = gridIndices[size_t(GridOffset(u, v))];
                const int i1 = gridIndices[size_t(GridOffset(u + 1, v))];
                const int i2 = gridIndices[size_t(GridOffset(u + 1, v + 1))];
                const int i3 = gridIndices[size_t(GridOffset(u, v + 1))];
                EmitQuadCellTriangles(mesh, i0, i1, i2, i3);
            }
        }
    }

    gFallbackSelectedMesh = nullptr;
    gFallbackCoarseFaceForPtex = nullptr;
    gFallbackFaceStartOffsets = nullptr;
    return mesh;
}

static bool ProcessSelectedMesh(const SelectedSubdivMesh &selected,
                                const UsdCameraInfo &usdCameraInfo,
                                const fs::path &outputDir,
                                EvaluationMode evaluationMode,
                                int uniformRate,
                                const char *farOutputName,
                                const char *nearOutputName)
{
    if (selected.faceVertexCounts.size() > 417)
    {
        const int faceVertexCount417 = selected.faceVertexCounts[417];
        printf("process_mesh_probe [%s] coarse_face=417 face_vertex_count=%d is_triangle=%d\n",
               selected.path.GetString().c_str(),
               faceVertexCount417,
               faceVertexCount417 == 3 ? 1 : 0);
    }

    Far::TopologyDescriptor desc = {};
    desc.numVertices = int(selected.points.size());
    desc.numFaces = int(selected.faceVertexCounts.size());
    desc.numVertsPerFace = selected.faceVertexCounts.cdata();
    desc.vertIndicesPerFace = selected.faceVertexIndices.cdata();
    std::vector<int> creaseVertexPairs;
    std::vector<float> creaseSharpnesses;
    BuildCreasePairs(selected, creaseVertexPairs, creaseSharpnesses);
    desc.numCreases = int(creaseSharpnesses.size());
    desc.creaseVertexIndexPairs = creaseVertexPairs.empty() ? nullptr : creaseVertexPairs.data();
    desc.creaseWeights = creaseSharpnesses.empty() ? nullptr : creaseSharpnesses.data();
    desc.numCorners = int(std::min(selected.cornerIndices.size(), selected.cornerSharpnesses.size()));
    desc.cornerVertexIndices = selected.cornerIndices.empty() ? nullptr : selected.cornerIndices.cdata();
    desc.cornerWeights = selected.cornerSharpnesses.empty() ? nullptr : selected.cornerSharpnesses.cdata();
    desc.numHoles = int(selected.holeIndices.size());
    desc.holeIndices = selected.holeIndices.empty() ? nullptr : selected.holeIndices.cdata();

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
    PrintRawPatchLookupProbe(patchMap, numPtexFaces, selected.path.GetString().c_str());

    std::vector<int> coarseFaceForPtex;
    std::vector<int> quadrantForPtex;
    BuildPtexMaps(selected, ptexIndices, numPtexFaces, coarseFaceForPtex, quadrantForPtex);
    std::vector<int> faceStartOffsets(selected.faceVertexCounts.size(), 0);
    int faceStart = 0;
    for (size_t face = 0; face < selected.faceVertexCounts.size(); face++)
    {
        faceStartOffsets[face] = faceStart;
        faceStart += selected.faceVertexCounts[face];
    }
    PrintNonEvaluableFaceInfo(
        patchMap, coarseFaceForPtex, selected, selected.path.GetString().c_str());
    PrintRawNonEvaluableFaceInfo(
        patchMap, coarseFaceForPtex, selected, selected.path.GetString().c_str());
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

    if (usdCameraInfo.found)
    {
        printf("Using USD camera: %s distance=%f\n",
               usdCameraInfo.path.GetString().c_str(),
               usdCameraInfo.distanceToMeshCenter);
    }
    else
    {
        printf("No USD camera metadata found in input; using fallback mesh-scaled camera distance.\n");
    }

    const std::array<CameraPreset, 2> cameras = BuildCameraPresets(selected, usdCameraInfo);
    const AdaptiveSettings settings = {};
    const char *outputNames[2] = {farOutputName, nearOutputName};
    const char *modeLabel =
        evaluationMode == EvaluationMode::Uniform ? "uniform_grid_no_stitch" : "adaptive_nosplit";
    for (int i = 0; i < 2; i++)
    {
        const CameraPreset &camera = cameras[size_t(i)];
        TessMesh tessMesh = {};
        if (evaluationMode == EvaluationMode::Uniform)
        {
            tessMesh = TessellateUniformNoStitch(patchMap,
                                                 *patchTable,
                                                 refinedPositions,
                                                 coarseFaceForPtex,
                                                 faceStartOffsets,
                                                 selected,
                                                 numPtexFaces,
                                                 uniformRate);
        }
        else
        {
            tessMesh = TessellateAdaptiveNoSplit(patchMap,
                                                 *patchTable,
                                                 refinedPositions,
                                                 adjacency,
                                                 coarseFaceForPtex,
                                                 faceStartOffsets,
                                                 selected,
                                                 numPtexFaces,
                                                 camera,
                                                 settings);
        }
        const fs::path outPath = outputDir / outputNames[i];
        const int uniformRateMetadata = evaluationMode == EvaluationMode::Uniform ? uniformRate : -1;
        if (!ybi::testio::WriteAdaptiveObj(tessMesh.positions,
                                           tessMesh.indices,
                                           outPath,
                                           selected,
                                           camera.name,
                                           camera.distanceToTarget,
                                           settings.numEdgeSamples,
                                           settings.pixelSpacing,
                                           settings.minRate,
                                           settings.maxRate,
                                           usdCameraInfo,
                                           modeLabel,
                                           uniformRateMetadata))
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

int main(int argc, char **argv)
{
    const CliOptions options = ParseCli(argc, argv);
    fs::create_directories(options.outputDir);

    SelectedSubdivMesh selected = {};
    UsdCameraInfo usdCameraInfo = {};
    if (!ybi::testio::LoadSelectedSubdivFromJson(options.inputJsonPath, selected, usdCameraInfo))
    {
        return 1;
    }

    printf("Loaded control cage JSON: %s\n", options.inputJsonPath.string().c_str());
    printf("  prim=%s scheme=%s vertices=%zu faces=%zu creases=%zu corners=%zu holes=%zu\n",
           selected.path.GetString().c_str(),
           selected.subdivisionScheme.c_str(),
           selected.points.size(),
           selected.faceVertexCounts.size(),
           selected.creaseSharpnesses.size(),
           selected.cornerIndices.size(),
           selected.holeIndices.size());

    if (!ProcessSelectedMesh(selected,
                             usdCameraInfo,
                             options.outputDir,
                             options.evaluationMode,
                             options.uniformRate,
                             options.farOutputName.c_str(),
                             options.nearOutputName.c_str()))
    {
        return 1;
    }

    return 0;
}
