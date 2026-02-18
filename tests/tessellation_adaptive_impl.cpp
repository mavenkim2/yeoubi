#include "bvh/usd_subdiv_select.h"
#include "bvh/usd_camera_utils.h"

#include <opensubdiv/far/patchMap.h>
#include <opensubdiv/far/patchTable.h>
#include <opensubdiv/far/patchTableFactory.h>
#include <opensubdiv/far/primvarRefiner.h>
#include <opensubdiv/far/ptexIndices.h>
#include <opensubdiv/far/topologyDescriptor.h>
#include <opensubdiv/sdc/types.h>

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <numeric>
#include <optional>
#include <sstream>
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

struct MeshValidationStats
{
    int outOfBoundsTriangles = 0;
    int degenerateIndexTriangles = 0;
    int degenerateAreaTriangles = 0;
    int maxVertexTriangleUse = 0;
    int suspiciousOverusedVertices = 0;
    int zeroPositionVertices = 0;
};

struct EvalDebugStats
{
    int patchLookupFailures = 0;
    int fallbackEvaluations = 0;
    int rawLookupMissesRecovered = 0;
};

static EvalDebugStats gEvalDebugStats = {};
static std::vector<std::string> gRecoveredRawMissExamples;
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

struct CliOptions
{
    fs::path inputJsonPath;
    fs::path outputDir = fs::path("tests") / "bvh" / "out";
    std::string controlCageOutputName = "subdiv_control_cage.obj";
    std::string farOutputName = "subdiv_limit_adaptive_nosplit_far.obj";
    std::string nearOutputName = "subdiv_limit_adaptive_nosplit_near.obj";
};

static std::string ExtractJsonString(const std::string &text, const std::string &key)
{
    const std::string token = "\"" + key + "\"";
    const size_t keyPos = text.find(token);
    if (keyPos == std::string::npos)
    {
        return {};
    }
    const size_t colonPos = text.find(':', keyPos);
    if (colonPos == std::string::npos)
    {
        return {};
    }
    const size_t quoteStart = text.find('"', colonPos + 1);
    if (quoteStart == std::string::npos)
    {
        return {};
    }
    const size_t quoteEnd = text.find('"', quoteStart + 1);
    if (quoteEnd == std::string::npos || quoteEnd <= quoteStart)
    {
        return {};
    }
    return text.substr(quoteStart + 1, quoteEnd - quoteStart - 1);
}

static std::optional<float> ExtractJsonFloat(const std::string &text, const std::string &key)
{
    const std::string token = "\"" + key + "\"";
    const size_t keyPos = text.find(token);
    if (keyPos == std::string::npos)
    {
        return std::nullopt;
    }
    const size_t colonPos = text.find(':', keyPos);
    if (colonPos == std::string::npos)
    {
        return std::nullopt;
    }

    size_t begin = colonPos + 1;
    while (begin < text.size() &&
           (text[begin] == ' ' || text[begin] == '\t' || text[begin] == '\n' || text[begin] == '\r'))
    {
        begin++;
    }
    if (begin >= text.size())
    {
        return std::nullopt;
    }

    size_t end = begin;
    while (end < text.size() &&
           (std::isdigit((unsigned char)text[end]) || text[end] == '.' || text[end] == '-' ||
            text[end] == '+' || text[end] == 'e' || text[end] == 'E'))
    {
        end++;
    }
    if (end <= begin)
    {
        return std::nullopt;
    }

    return std::stof(text.substr(begin, end - begin));
}

static std::string ExtractJsonArray(const std::string &text, const std::string &key)
{
    const std::string token = "\"" + key + "\"";
    const size_t keyPos = text.find(token);
    if (keyPos == std::string::npos)
    {
        return {};
    }
    const size_t bracketStart = text.find('[', keyPos);
    if (bracketStart == std::string::npos)
    {
        return {};
    }

    int depth = 0;
    for (size_t i = bracketStart; i < text.size(); i++)
    {
        if (text[i] == '[')
        {
            depth++;
        }
        else if (text[i] == ']')
        {
            depth--;
            if (depth == 0)
            {
                return text.substr(bracketStart, i - bracketStart + 1);
            }
        }
    }

    return {};
}

static std::vector<float> ParseFloatArray(const std::string &arrayText)
{
    std::vector<float> values;
    size_t index = 0;
    while (index < arrayText.size())
    {
        while (index < arrayText.size() &&
               !(std::isdigit((unsigned char)arrayText[index]) || arrayText[index] == '-' ||
                 arrayText[index] == '+' || arrayText[index] == '.'))
        {
            index++;
        }
        if (index >= arrayText.size())
        {
            break;
        }

        size_t endIndex = index + 1;
        while (endIndex < arrayText.size() &&
               (std::isdigit((unsigned char)arrayText[endIndex]) || arrayText[endIndex] == '.' ||
                arrayText[endIndex] == '-' || arrayText[endIndex] == '+' ||
                arrayText[endIndex] == 'e' || arrayText[endIndex] == 'E'))
        {
            endIndex++;
        }

        values.push_back(std::stof(arrayText.substr(index, endIndex - index)));
        index = endIndex;
    }
    return values;
}

static std::vector<int> ParseIntArray(const std::string &arrayText)
{
    const std::vector<float> parsed = ParseFloatArray(arrayText);
    std::vector<int> values;
    values.reserve(parsed.size());
    for (size_t i = 0; i < parsed.size(); i++)
    {
        values.push_back((int)parsed[i]);
    }
    return values;
}

static void PrintUsage(const char *exeName)
{
    printf("Usage: %s --input-json <path> [--out-dir <dir>] [--control-cage-out <obj>] "
           "[--far-out <obj>] [--near-out <obj>]\n",
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
        if (arg == "--control-cage-out")
        {
            if (i + 1 >= argc)
            {
                PrintUsage(argv[0]);
                std::abort();
            }
            options.controlCageOutputName = argv[++i];
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

static bool LoadSelectedSubdivFromJson(const fs::path &jsonPath,
                                       SelectedSubdivMesh &selectedOut,
                                       UsdCameraInfo &usdCameraInfoOut)
{
    std::ifstream input(jsonPath, std::ios::in | std::ios::binary);
    if (!input.is_open())
    {
        fprintf(stderr, "Failed to open input JSON: %s\n", jsonPath.string().c_str());
        return false;
    }
    const std::string json((std::istreambuf_iterator<char>(input)), std::istreambuf_iterator<char>());

    const std::string sourcePrim = ExtractJsonString(json, "source_prim");
    if (!sourcePrim.empty())
    {
        selectedOut.path = pxr::SdfPath(sourcePrim);
    }
    selectedOut.subdivisionScheme = ExtractJsonString(json, "scheme");
    if (selectedOut.subdivisionScheme.empty())
    {
        selectedOut.subdivisionScheme = "catmullClark";
    }

    const std::vector<float> pointScalars = ParseFloatArray(ExtractJsonArray(json, "points"));
    if (pointScalars.size() % 3 != 0)
    {
        fprintf(stderr, "Invalid points array in JSON (not xyz triplets): %s\n", jsonPath.string().c_str());
        return false;
    }
    selectedOut.points.resize(pointScalars.size() / 3);
    for (size_t i = 0; i + 2 < pointScalars.size(); i += 3)
    {
        selectedOut.points[i / 3] = pxr::GfVec3f(pointScalars[i], pointScalars[i + 1], pointScalars[i + 2]);
    }

    const std::vector<int> faceVertexCounts = ParseIntArray(ExtractJsonArray(json, "face_vertex_counts"));
    const std::vector<int> faceVertexIndices = ParseIntArray(ExtractJsonArray(json, "face_vertex_indices"));
    selectedOut.faceVertexCounts = pxr::VtIntArray(faceVertexCounts.begin(), faceVertexCounts.end());
    selectedOut.faceVertexIndices = pxr::VtIntArray(faceVertexIndices.begin(), faceVertexIndices.end());

    const std::vector<int> cornerIndices = ParseIntArray(ExtractJsonArray(json, "corner_indices"));
    const std::vector<float> cornerSharpnesses =
        ParseFloatArray(ExtractJsonArray(json, "corner_sharpnesses"));
    const std::vector<int> creaseIndices = ParseIntArray(ExtractJsonArray(json, "crease_indices"));
    const std::vector<int> creaseLengths = ParseIntArray(ExtractJsonArray(json, "crease_lengths"));
    const std::vector<float> creaseSharpnesses =
        ParseFloatArray(ExtractJsonArray(json, "crease_sharpnesses"));
    const std::vector<int> holeIndices = ParseIntArray(ExtractJsonArray(json, "hole_indices"));

    selectedOut.cornerIndices = pxr::VtIntArray(cornerIndices.begin(), cornerIndices.end());
    selectedOut.cornerSharpnesses =
        pxr::VtFloatArray(cornerSharpnesses.begin(), cornerSharpnesses.end());
    selectedOut.creaseIndices = pxr::VtIntArray(creaseIndices.begin(), creaseIndices.end());
    selectedOut.creaseLengths = pxr::VtIntArray(creaseLengths.begin(), creaseLengths.end());
    selectedOut.creaseSharpnesses =
        pxr::VtFloatArray(creaseSharpnesses.begin(), creaseSharpnesses.end());
    selectedOut.holeIndices = pxr::VtIntArray(holeIndices.begin(), holeIndices.end());

    const int expectedFaceIndexCount =
        std::accumulate(selectedOut.faceVertexCounts.begin(), selectedOut.faceVertexCounts.end(), 0);
    if (expectedFaceIndexCount != int(selectedOut.faceVertexIndices.size()))
    {
        fprintf(stderr,
                "Invalid face topology in JSON: expected %d face indices, found %zu (%s)\n",
                expectedFaceIndexCount,
                selectedOut.faceVertexIndices.size(),
                jsonPath.string().c_str());
        return false;
    }

    usdCameraInfoOut = {};
    const std::string cameraPath = ExtractJsonString(json, "usd_camera_path");
    const std::optional<float> cameraDistance =
        ExtractJsonFloat(json, "usd_camera_distance_to_mesh_center");
    if (!cameraPath.empty() && cameraDistance.has_value())
    {
        usdCameraInfoOut.found = true;
        usdCameraInfoOut.path = pxr::SdfPath(cameraPath);
        usdCameraInfoOut.distanceToMeshCenter = cameraDistance.value();
    }

    return true;
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
    if (!rawHandle && handle)
    {
        gEvalDebugStats.rawLookupMissesRecovered++;
        if (gRecoveredRawMissExamples.size() < 16)
        {
            char line[256];
            snprintf(line,
                     sizeof(line),
                     "ptex=%d req=(%.6f,%.6f) recovered=(%.6f,%.6f)",
                     ptexFace,
                     u,
                     v,
                     evalU,
                     evalV);
            gRecoveredRawMissExamples.push_back(std::string(line));
        }
    }
    if (!handle)
    {
        gEvalDebugStats.patchLookupFailures++;
        const pxr::GfVec3f coarseFallback = EvaluateCoarseFallbackPosition(ptexFace, u, v);
        if (!(fabsf(coarseFallback[0]) < 1e-20f && fabsf(coarseFallback[1]) < 1e-20f &&
              fabsf(coarseFallback[2]) < 1e-20f))
        {
            gEvalDebugStats.fallbackEvaluations++;
            return coarseFallback;
        }
        if (!positions.empty())
        {
            gEvalDebugStats.fallbackEvaluations++;
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

static MeshValidationStats ValidateMesh(const TessMesh &mesh, const char *label)
{
    MeshValidationStats stats = {};
    if (mesh.positions.empty() || mesh.indices.empty())
    {
        return stats;
    }

    std::vector<int> usage(mesh.positions.size(), 0);
    for (const pxr::GfVec3f &p : mesh.positions)
    {
        if (fabsf(p[0]) < 1e-20f && fabsf(p[1]) < 1e-20f && fabsf(p[2]) < 1e-20f)
        {
            stats.zeroPositionVertices++;
        }
    }
    for (size_t i = 0; i + 2 < mesh.indices.size(); i += 3)
    {
        const int i0 = mesh.indices[i + 0];
        const int i1 = mesh.indices[i + 1];
        const int i2 = mesh.indices[i + 2];
        const bool outOfBounds = i0 < 0 || i1 < 0 || i2 < 0 || i0 >= int(mesh.positions.size()) ||
                                 i1 >= int(mesh.positions.size()) || i2 >= int(mesh.positions.size());
        if (outOfBounds)
        {
            stats.outOfBoundsTriangles++;
            continue;
        }

        if (i0 == i1 || i1 == i2 || i2 == i0)
        {
            stats.degenerateIndexTriangles++;
            continue;
        }

        const float areaSq = TriangleAreaSquared(
            mesh.positions[size_t(i0)], mesh.positions[size_t(i1)], mesh.positions[size_t(i2)]);
        if (areaSq < 1e-18f)
        {
            stats.degenerateAreaTriangles++;
            continue;
        }

        usage[size_t(i0)]++;
        usage[size_t(i1)]++;
        usage[size_t(i2)]++;
    }

    long long totalUse = 0;
    for (int count : usage)
    {
        stats.maxVertexTriangleUse = std::max(stats.maxVertexTriangleUse, count);
        totalUse += count;
    }
    const float meanUse = usage.empty() ? 0.0f : float(totalUse) / float(usage.size());
    for (int count : usage)
    {
        if (count > 500 && count > int(meanUse * 50.0f))
        {
            stats.suspiciousOverusedVertices++;
        }
    }

    printf("Mesh validation [%s]: oob=%d deg_idx=%d deg_area=%d max_use=%d suspicious_use=%d zero_pos=%d\n",
           label,
           stats.outOfBoundsTriangles,
           stats.degenerateIndexTriangles,
           stats.degenerateAreaTriangles,
           stats.maxVertexTriangleUse,
           stats.suspiciousOverusedVertices,
           stats.zeroPositionVertices);
    return stats;
}

static bool CleanupInvalidTriangles(TessMesh &mesh)
{
    bool changed = false;
    std::vector<int> cleaned;
    cleaned.reserve(mesh.indices.size());

    for (size_t i = 0; i + 2 < mesh.indices.size(); i += 3)
    {
        const int i0 = mesh.indices[i + 0];
        const int i1 = mesh.indices[i + 1];
        const int i2 = mesh.indices[i + 2];
        const bool outOfBounds = i0 < 0 || i1 < 0 || i2 < 0 || i0 >= int(mesh.positions.size()) ||
                                 i1 >= int(mesh.positions.size()) || i2 >= int(mesh.positions.size());
        if (outOfBounds || IsTriangleDegenerate(mesh, i0, i1, i2))
        {
            changed = true;
            continue;
        }

        cleaned.push_back(i0);
        cleaned.push_back(i1);
        cleaned.push_back(i2);
    }

    if (changed)
    {
        mesh.indices.swap(cleaned);
    }
    return changed;
}

static int FixZeroPositionVertices(TessMesh &mesh)
{
    if (mesh.positions.empty() || mesh.indices.empty())
    {
        return 0;
    }

    std::vector<pxr::GfVec3f> sum(mesh.positions.size(), pxr::GfVec3f(0.0f));
    std::vector<int> count(mesh.positions.size(), 0);

    for (size_t i = 0; i + 2 < mesh.indices.size(); i += 3)
    {
        const int idx[3] = {mesh.indices[i + 0], mesh.indices[i + 1], mesh.indices[i + 2]};
        for (int a = 0; a < 3; a++)
        {
            if (idx[a] < 0 || idx[a] >= int(mesh.positions.size()))
            {
                continue;
            }
            for (int b = 0; b < 3; b++)
            {
                if (a == b || idx[b] < 0 || idx[b] >= int(mesh.positions.size()))
                {
                    continue;
                }
                sum[size_t(idx[a])] += mesh.positions[size_t(idx[b])];
                count[size_t(idx[a])]++;
            }
        }
    }

    pxr::GfVec3f fallback(0.0f);
    for (const pxr::GfVec3f &p : mesh.positions)
    {
        if (!(fabsf(p[0]) < 1e-20f && fabsf(p[1]) < 1e-20f && fabsf(p[2]) < 1e-20f))
        {
            fallback = p;
            break;
        }
    }

    int fixed = 0;
    for (size_t i = 0; i < mesh.positions.size(); i++)
    {
        const pxr::GfVec3f p = mesh.positions[i];
        const bool isZero = fabsf(p[0]) < 1e-20f && fabsf(p[1]) < 1e-20f && fabsf(p[2]) < 1e-20f;
        if (!isZero)
        {
            continue;
        }

        if (count[i] > 0)
        {
            mesh.positions[i] = sum[i] / float(count[i]);
            fixed++;
        }
        else
        {
            mesh.positions[i] = fallback;
            fixed++;
        }
    }
    return fixed;
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
    for (int face = 0; face < numPtexFaces; face++)
    {
        faceIsEvaluable[size_t(face)] = FaceHasPatchCoverage(patchMap, face);
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
                if (!IsTriangleDegenerate(mesh, i0, i1, i2))
                {
                    mesh.indices.push_back(i0);
                    mesh.indices.push_back(i1);
                    mesh.indices.push_back(i2);
                }
                if (!IsTriangleDegenerate(mesh, i0, i2, i3))
                {
                    mesh.indices.push_back(i0);
                    mesh.indices.push_back(i2);
                    mesh.indices.push_back(i3);
                }
            }
        }
    }

    gFallbackSelectedMesh = nullptr;
    gFallbackCoarseFaceForPtex = nullptr;
    gFallbackFaceStartOffsets = nullptr;

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
        out << "f";
        for (int i = 0; i < faceVertexCount; i++)
        {
            const int index = source.faceVertexIndices[indexOffset + i];
            out << " " << (index + 1);
        }
        out << "\n";
        indexOffset += faceVertexCount;
    }

    return out.good();
}

static bool ProcessSelectedMesh(const SelectedSubdivMesh &selected,
                                const UsdCameraInfo &usdCameraInfo,
                                const fs::path &outputDir,
                                const char *controlCageName,
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
    // Temporary debug mode: disable crease/corner/hole handling to isolate artifacts.
    desc.numCreases = 0;
    desc.creaseVertexIndexPairs = nullptr;
    desc.creaseWeights = nullptr;
    desc.numCorners = 0;
    desc.cornerVertexIndices = nullptr;
    desc.cornerWeights = nullptr;
    desc.numHoles = 0;
    desc.holeIndices = nullptr;

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
    printf("Control-cage features: crease/corner/hole handling disabled for tessellation debug.\n");

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
    for (int i = 0; i < 2; i++)
    {
        const CameraPreset &camera = cameras[size_t(i)];
        gEvalDebugStats = {};
        gRecoveredRawMissExamples.clear();
        TessMesh tessMesh = TessellateAdaptiveNoSplit(
            patchMap,
            *patchTable,
            refinedPositions,
            adjacency,
            coarseFaceForPtex,
            faceStartOffsets,
            selected,
            numPtexFaces,
            camera,
            settings);
        for (int pass = 0; pass < 64; pass++)
        {
            const MeshValidationStats stats =
                ValidateMesh(tessMesh, camera.name);
            const bool hasCriticalInvalid = stats.outOfBoundsTriangles > 0 ||
                                            stats.degenerateIndexTriangles > 0 ||
                                            stats.degenerateAreaTriangles > 0;
            const bool hasZeroPositions = stats.zeroPositionVertices > 0;
            if (!hasCriticalInvalid && !hasZeroPositions)
            {
                break;
            }

            bool changed = false;
            if (hasCriticalInvalid)
            {
                changed = CleanupInvalidTriangles(tessMesh) || changed;
            }
            if (hasZeroPositions)
            {
                const int fixedZeroCount = FixZeroPositionVertices(tessMesh);
                changed = fixedZeroCount > 0 || changed;
                printf("Zero-fix [%s]: fixed=%d pass=%d\n", camera.name, fixedZeroCount, pass);
            }
            if (!changed)
            {
                break;
            }
        }
        printf("Eval debug [%s]: patch_lookup_failures=%d fallbacks=%d\n",
               camera.name,
               gEvalDebugStats.patchLookupFailures,
               gEvalDebugStats.fallbackEvaluations);
        printf("Eval debug [%s]: raw_misses_recovered=%d\n",
               camera.name,
               gEvalDebugStats.rawLookupMissesRecovered);
        for (const std::string &example : gRecoveredRawMissExamples)
        {
            printf("  recovered_raw_miss [%s] %s\n", camera.name, example.c_str());
        }
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

int main(int argc, char **argv)
{
    const CliOptions options = ParseCli(argc, argv);
    fs::create_directories(options.outputDir);

    SelectedSubdivMesh selected = {};
    UsdCameraInfo usdCameraInfo = {};
    if (!LoadSelectedSubdivFromJson(options.inputJsonPath, selected, usdCameraInfo))
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
                             options.controlCageOutputName.c_str(),
                             options.farOutputName.c_str(),
                             options.nearOutputName.c_str()))
    {
        return 1;
    }

    return 0;
}
