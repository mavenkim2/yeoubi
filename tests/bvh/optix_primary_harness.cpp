#include "device/cuda_device.h"
#include "io/usd/load.h"
#include "scene/scene.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "third_party/stb_image_write.h"
#include "util/array.h"
#include "util/float3.h"
#include "util/float3x4.h"
#include "util/float4.h"
#include "util/float4x4.h"
#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <optional>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include <optix_stubs.h>

using namespace ybi;

#ifndef YBI_OPTIX_PRIMARY_PTX_PATH
#define YBI_OPTIX_PRIMARY_PTX_PATH ""
#endif

namespace
{
static void OptixCheckImpl(OptixResult result, const char *expr, const char *file, int line)
{
    if (result != OPTIX_SUCCESS)
    {
        fprintf(stderr,
                "OptiX call failed: %s -> %s (%s) at %s:%d\n",
                expr,
                optixGetErrorName(result),
                optixGetErrorString(result),
                file,
                line);
        std::abort();
    }
}

#define OPTIX_CHECK(expr) OptixCheckImpl((expr), #expr, __FILE__, __LINE__)

struct LaunchParams
{
    struct InstanceGeomRef
    {
        unsigned long long positions;
        unsigned long long indices;
        int numPositions;
        int numIndices;
    };

    struct WireframeConfig
    {
        float lineWidth;
        float lineFeather;
        float edgeDarkness;
        float padding;
    };

    OptixTraversableHandle traversable;
    CUdeviceptr image;
    int width;
    int height;
    ybi::float3 cameraOrigin;
    ybi::float3 cameraU;
    ybi::float3 cameraV;
    ybi::float3 cameraW;
    WireframeConfig wireframe;
    int integrator;
    int spp;
    float aoBias;
    float aoMaxDistance;
    unsigned long long instanceGeomRefs;
    int instanceGeomRefCount;
};

enum class RenderType
{
    Triangle,
    Cluster,
    Curve,
    USD
};

enum class IntegratorType
{
    Primary,
    AO
};

struct CliOptions
{
    RenderType type = RenderType::Triangle;
    IntegratorType integrator = IntegratorType::Primary;
    std::string inputPath = "tests/bvh/out/stoat_body_selected.obj";
    std::string outputPath = "optix_triangle_gas.png";
    std::optional<ybi::float3> cameraPosition;
    std::optional<ybi::float3> lookAt;
    int spp = 1;
};

struct RenderCameraOverride
{
    int width = 1280;
    int height = 720;
    ybi::float3 origin = ybi::make_float3(0.0f, 0.0f, 0.0f);
    ybi::float3 U = ybi::make_float3(1.0f, 0.0f, 0.0f);
    ybi::float3 V = ybi::make_float3(0.0f, 1.0f, 0.0f);
    ybi::float3 W = ybi::make_float3(0.0f, 0.0f, 1.0f);
};

static std::string ReadTextFile(const std::string &path)
{
    std::ifstream input(path, std::ios::in | std::ios::binary);
    if (!input.is_open())
    {
        fprintf(stderr, "Failed to open PTX file: %s\n", path.c_str());
        std::abort();
    }
    return std::string((std::istreambuf_iterator<char>(input)), std::istreambuf_iterator<char>());
}

static ybi::float3 Cross(const ybi::float3 &a, const ybi::float3 &b)
{
    return ybi::make_float3(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
}

static ybi::float3 Normalize(const ybi::float3 &v)
{
    const float lenSq = ybi::dot(v, v);
    if (lenSq <= 0.0f)
    {
        return ybi::make_float3(0.0f, 0.0f, 1.0f);
    }
    const float invLen = 1.0f / std::sqrt(lenSq);
    return v * invLen;
}

static bool InvertAffine(const ybi::float4x4 &m, ybi::float4x4 &out)
{
    const float a00 = m.m[0][0], a01 = m.m[0][1], a02 = m.m[0][2];
    const float a10 = m.m[1][0], a11 = m.m[1][1], a12 = m.m[1][2];
    const float a20 = m.m[2][0], a21 = m.m[2][1], a22 = m.m[2][2];

    const float c00 = a11 * a22 - a12 * a21;
    const float c01 = a02 * a21 - a01 * a22;
    const float c02 = a01 * a12 - a02 * a11;
    const float c10 = a12 * a20 - a10 * a22;
    const float c11 = a00 * a22 - a02 * a20;
    const float c12 = a02 * a10 - a00 * a12;
    const float c20 = a10 * a21 - a11 * a20;
    const float c21 = a01 * a20 - a00 * a21;
    const float c22 = a00 * a11 - a01 * a10;

    const float det = a00 * c00 + a01 * c10 + a02 * c20;
    if (std::abs(det) < 1e-8f)
    {
        return false;
    }
    const float invDet = 1.0f / det;

    const float r00 = c00 * invDet, r01 = c01 * invDet, r02 = c02 * invDet;
    const float r10 = c10 * invDet, r11 = c11 * invDet, r12 = c12 * invDet;
    const float r20 = c20 * invDet, r21 = c21 * invDet, r22 = c22 * invDet;

    const float tx = m.m[0][3];
    const float ty = m.m[1][3];
    const float tz = m.m[2][3];
    const float itx = -(r00 * tx + r01 * ty + r02 * tz);
    const float ity = -(r10 * tx + r11 * ty + r12 * tz);
    const float itz = -(r20 * tx + r21 * ty + r22 * tz);

    out = ybi::float4x4(
        r00, r01, r02, itx, r10, r11, r12, ity, r20, r21, r22, itz, 0.0f, 0.0f, 0.0f, 1.0f);
    return true;
}

static std::optional<RenderCameraOverride> BuildUsdRenderCamera(const Camera &camera)
{
    if (camera.viewportWidth <= 0 || camera.viewportHeight <= 0)
    {
        return std::nullopt;
    }

    ybi::float4x4 worldFromCamera = {};
    if (!InvertAffine(camera.cameraFromWorld, worldFromCamera))
    {
        return std::nullopt;
    }

    ybi::float3 right = Normalize(ybi::make_float3(
        worldFromCamera.m[0][0], worldFromCamera.m[1][0], worldFromCamera.m[2][0]));
    ybi::float3 up = Normalize(ybi::make_float3(
        worldFromCamera.m[0][1], worldFromCamera.m[1][1], worldFromCamera.m[2][1]));
    ybi::float3 forward = Normalize(ybi::make_float3(
        -worldFromCamera.m[0][2], -worldFromCamera.m[1][2], -worldFromCamera.m[2][2]));
    if (std::abs(ybi::dot(forward, forward)) < 1e-6f)
    {
        return std::nullopt;
    }

    const float m00 = camera.clipFromCamera.m[0][0];
    const float m11 = camera.clipFromCamera.m[1][1];
    float tanHalfFov = std::tan(45.0f * 0.5f * 3.14159265358979323846f / 180.0f);
    float aspect =
        static_cast<float>(camera.viewportWidth) / static_cast<float>(camera.viewportHeight);
    if (std::abs(m00) > 1e-6f && std::abs(m11) > 1e-6f)
    {
        tanHalfFov = 1.0f / std::abs(m11);
        aspect = std::abs(m11 / m00);
    }
    if (m11 < 0.0f)
    {
        up = up * -1.0f;
    }

    RenderCameraOverride out = {};
    out.width = camera.viewportWidth;
    out.height = camera.viewportHeight;
    out.origin = ybi::make_float3(
        worldFromCamera.m[0][3], worldFromCamera.m[1][3], worldFromCamera.m[2][3]);
    out.U = right * (aspect * tanHalfFov);
    out.V = up * tanHalfFov;
    out.W = forward;
    return out;
}

static bool ParseFloat3(int argc, char **argv, int startIndex, ybi::float3 &valueOut)
{
    if (startIndex + 2 >= argc)
    {
        return false;
    }
    valueOut.x = std::stof(argv[startIndex + 0]);
    valueOut.y = std::stof(argv[startIndex + 1]);
    valueOut.z = std::stof(argv[startIndex + 2]);
    return true;
}

static void PrintUsage(const char *exeName)
{
    printf("Usage: %s [--type triangle|cluster|curve|usd] [--file path] [--out path] "
           "[--integrator primary|ao] [--spp N] [--cam-pos x y z] [--look-at x y z]\n",
           exeName);
    printf("  --type triangle|cluster|curve|usd\n");
    printf("  --file OBJ path for triangle/cluster; JSON path for curve; USDA/USD path for usd\n");
    printf("  --out PNG output path\n");
    printf("  --integrator primary|ao\n");
    printf("  --spp samples per pixel for ao integrator\n");
    printf("  --cam-pos optional camera position override\n");
    printf("  --look-at optional look-at target (default bounds center)\n");
}

static CliOptions ParseCli(int argc, char **argv)
{
    CliOptions options = {};
    for (int i = 1; i < argc; i++)
    {
        const std::string arg = argv[i];
        if (arg == "--type")
        {
            if (i + 1 >= argc)
            {
                PrintUsage(argv[0]);
                std::abort();
            }
            const std::string value = argv[++i];
            if (value == "triangle")
            {
                options.type = RenderType::Triangle;
                if (options.outputPath == "optix_triangle_gas.png")
                {
                    options.outputPath = "optix_triangle_gas.png";
                }
            }
            else if (value == "cluster")
            {
                options.type = RenderType::Cluster;
                options.outputPath = "optix_cluster_gas.png";
            }
            else if (value == "curve")
            {
                options.type = RenderType::Curve;
                options.inputPath = "tests/bvh/out/selected_curve.json";
                options.outputPath = "optix_curve_gas.png";
            }
            else if (value == "usd")
            {
                options.type = RenderType::USD;
                options.outputPath = "optix_usd_scene.png";
            }
            else
            {
                PrintUsage(argv[0]);
                std::abort();
            }
            continue;
        }
        if (arg == "--file")
        {
            if (i + 1 >= argc)
            {
                PrintUsage(argv[0]);
                std::abort();
            }
            options.inputPath = argv[++i];
            continue;
        }
        if (arg == "--integrator")
        {
            if (i + 1 >= argc)
            {
                PrintUsage(argv[0]);
                std::abort();
            }
            const std::string value = argv[++i];
            if (value == "primary")
            {
                options.integrator = IntegratorType::Primary;
            }
            else if (value == "ao")
            {
                options.integrator = IntegratorType::AO;
            }
            else
            {
                PrintUsage(argv[0]);
                std::abort();
            }
            continue;
        }
        if (arg == "--spp")
        {
            if (i + 1 >= argc)
            {
                PrintUsage(argv[0]);
                std::abort();
            }
            options.spp = std::max(1, std::stoi(argv[++i]));
            continue;
        }
        if (arg == "--out")
        {
            if (i + 1 >= argc)
            {
                PrintUsage(argv[0]);
                std::abort();
            }
            options.outputPath = argv[++i];
            continue;
        }
        if (arg == "--cam-pos")
        {
            ybi::float3 value = {};
            if (!ParseFloat3(argc, argv, i + 1, value))
            {
                PrintUsage(argv[0]);
                std::abort();
            }
            options.cameraPosition = value;
            i += 3;
            continue;
        }
        if (arg == "--look-at")
        {
            ybi::float3 value = {};
            if (!ParseFloat3(argc, argv, i + 1, value))
            {
                PrintUsage(argv[0]);
                std::abort();
            }
            options.lookAt = value;
            i += 3;
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

    if (options.integrator == IntegratorType::Primary)
    {
        options.spp = 1;
    }

    return options;
}

static bool SavePNG(const char *filePath, const std::vector<uint8_t> &rgba, int width, int height)
{
    const int strideInBytes = width * 4;
    return stbi_write_png(filePath, width, height, 4, rgba.data(), strideInBytes) != 0;
}

static bool ParseObjIndex(const std::string &token, int &indexOut)
{
    const size_t slashPos = token.find('/');
    const std::string indexText =
        slashPos == std::string::npos ? token : token.substr(0, slashPos);
    if (indexText.empty())
    {
        return false;
    }
    indexOut = std::stoi(indexText);
    return true;
}

static Mesh LoadObjMesh(const std::string &path)
{
    std::ifstream input(path);
    if (!input.is_open())
    {
        fprintf(stderr, "Failed to open OBJ: %s\n", path.c_str());
        std::abort();
    }

    std::vector<ybi::float3> positions;
    std::vector<int> indices;
    std::string line;
    while (std::getline(input, line))
    {
        if (line.empty() || line[0] == '#')
        {
            continue;
        }
        std::istringstream ss(line);
        std::string tag;
        ss >> tag;
        if (tag == "v")
        {
            float x, y, z;
            ss >> x >> y >> z;
            positions.push_back(ybi::make_float3(x, y, z));
            continue;
        }
        if (tag == "f")
        {
            std::vector<int> face;
            std::string token;
            while (ss >> token)
            {
                int idx = 0;
                if (!ParseObjIndex(token, idx))
                {
                    continue;
                }
                if (idx < 0)
                {
                    idx = static_cast<int>(positions.size()) + idx;
                }
                else
                {
                    idx = idx - 1;
                }
                face.push_back(idx);
            }

            if (face.size() >= 3)
            {
                const int i0 = face[0];
                for (size_t i = 1; i + 1 < face.size(); i++)
                {
                    indices.push_back(i0);
                    indices.push_back(face[i]);
                    indices.push_back(face[i + 1]);
                }
            }
        }
    }

    Array<ybi::float3> meshPositions(positions);
    Array<int> meshIndices(indices);
    return Mesh(std::move(meshPositions), std::move(meshIndices));
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

static Curves LoadCurveJson(const std::string &path)
{
    std::ifstream input(path, std::ios::in | std::ios::binary);
    if (!input.is_open())
    {
        fprintf(stderr, "Failed to open curve JSON: %s\n", path.c_str());
        return {};
    }

    const std::string json((std::istreambuf_iterator<char>(input)),
                           std::istreambuf_iterator<char>());
    const std::string pointsArray = ExtractJsonArray(json, "points");
    const std::string widthsArray = ExtractJsonArray(json, "widths");
    const std::string curveCountsArray = ExtractJsonArray(json, "curve_vertex_counts");
    if (pointsArray.empty() || curveCountsArray.empty())
    {
        fprintf(stderr, "Curve JSON missing required arrays: %s\n", path.c_str());
        return {};
    }

    const std::vector<float> pointScalars = ParseFloatArray(pointsArray);
    if (pointScalars.size() % 3 != 0)
    {
        fprintf(stderr, "Curve JSON points are not xyz triplets: %s\n", path.c_str());
        return {};
    }

    std::vector<ybi::float3> positions;
    positions.reserve(pointScalars.size() / 3);
    for (size_t i = 0; i + 2 < pointScalars.size(); i += 3)
    {
        positions.push_back(
            ybi::make_float3(pointScalars[i + 0], pointScalars[i + 1], pointScalars[i + 2]));
    }

    const std::vector<int> curveVertexCounts = ParseIntArray(curveCountsArray);
    int expectedVertices = 0;
    for (size_t i = 0; i < curveVertexCounts.size(); i++)
    {
        expectedVertices += curveVertexCounts[i];
    }
    if ((size_t)expectedVertices != positions.size())
    {
        fprintf(stderr,
                "Curve JSON vertex count mismatch: counts=%d points=%zu in %s\n",
                expectedVertices,
                positions.size(),
                path.c_str());
        return {};
    }

    std::vector<float> widths = ParseFloatArray(widthsArray);
    if (widths.size() == 0)
    {
        widths.resize(positions.size(), 0.01f);
    }
    else if (widths.size() == 1)
    {
        widths.resize(positions.size(), widths[0]);
    }
    else if (widths.size() == curveVertexCounts.size())
    {
        std::vector<float> expanded;
        expanded.reserve(positions.size());
        for (size_t curveIndex = 0; curveIndex < curveVertexCounts.size(); curveIndex++)
        {
            for (int i = 0; i < curveVertexCounts[curveIndex]; i++)
            {
                expanded.push_back(widths[curveIndex]);
            }
        }
        widths = std::move(expanded);
    }
    else if (widths.size() != positions.size())
    {
        fprintf(stderr, "Unsupported widths count in curve JSON: %s\n", path.c_str());
        return {};
    }

    // USD-exported fur curve widths are often too tiny for this debug camera/resolution.
    // Normalize to a visible minimum and apply a small scale boost for preview rendering.
    const float widthScale = 1.0f;
    const float minWidth = 0.0f;
    for (size_t i = 0; i < widths.size(); i++)
    {
        widths[i] = std::max(widths[i] * widthScale, minWidth);
    }
    float widthMin = std::numeric_limits<float>::max();
    float widthMax = 0.0f;
    for (size_t i = 0; i < widths.size(); i++)
    {
        widthMin = std::min(widthMin, widths[i]);
        widthMax = std::max(widthMax, widths[i]);
    }
    printf("Curve JSON parsed: curves=%zu vertices=%zu widths=%zu widthMin=%f widthMax=%f\n",
           curveVertexCounts.size(),
           positions.size(),
           widths.size(),
           widthMin,
           widthMax);

    std::vector<int> curveVertexOffsets;
    curveVertexOffsets.reserve(curveVertexCounts.size());
    int offset = 0;
    for (size_t i = 0; i < curveVertexCounts.size(); i++)
    {
        curveVertexOffsets.push_back(offset);
        offset += curveVertexCounts[i];
    }

    Array<ybi::float3> curvePositions(positions);
    Array<float> curveWidths(widths);
    Array<int> offsets(curveVertexOffsets);
    return Curves(std::move(curvePositions), std::move(curveWidths), std::move(offsets));
}

static void ComputeBounds(const Mesh &mesh, ybi::float3 &boundsMinOut, ybi::float3 &boundsMaxOut)
{
    ybi::float3 boundsMin = ybi::make_float3(std::numeric_limits<float>::max());
    ybi::float3 boundsMax = ybi::make_float3(-std::numeric_limits<float>::max());
    for (size_t i = 0; i < mesh.positions.size(); i++)
    {
        const ybi::float3 p = mesh.positions[i];
        boundsMin.x = std::min(boundsMin.x, p.x);
        boundsMin.y = std::min(boundsMin.y, p.y);
        boundsMin.z = std::min(boundsMin.z, p.z);
        boundsMax.x = std::max(boundsMax.x, p.x);
        boundsMax.y = std::max(boundsMax.y, p.y);
        boundsMax.z = std::max(boundsMax.z, p.z);
    }
    boundsMinOut = boundsMin;
    boundsMaxOut = boundsMax;
}

static void
ComputeBounds(const Curves &curves, ybi::float3 &boundsMinOut, ybi::float3 &boundsMaxOut)
{
    ybi::float3 boundsMin = ybi::make_float3(std::numeric_limits<float>::max());
    ybi::float3 boundsMax = ybi::make_float3(-std::numeric_limits<float>::max());
    const Array<ybi::float3> &positions = curves.GetVertices();
    for (size_t i = 0; i < positions.size(); i++)
    {
        const ybi::float3 p = positions[i];
        boundsMin.x = std::min(boundsMin.x, p.x);
        boundsMin.y = std::min(boundsMin.y, p.y);
        boundsMin.z = std::min(boundsMin.z, p.z);
        boundsMax.x = std::max(boundsMax.x, p.x);
        boundsMax.y = std::max(boundsMax.y, p.y);
        boundsMax.z = std::max(boundsMax.z, p.z);
    }
    boundsMinOut = boundsMin;
    boundsMaxOut = boundsMax;
}

struct MeshAccelData
{
    OptixTraversableHandle gasHandle = {};
    CUdeviceptr positionsBuffer = 0;
    CUdeviceptr indicesBuffer = 0;
    int numPositions = 0;
    int numIndices = 0;
    ybi::float3 boundsMin = {};
    ybi::float3 boundsMax = {};
};

struct CurveAccelData
{
    OptixTraversableHandle gasHandle = {};
    ybi::float3 boundsMin = {};
    ybi::float3 boundsMax = {};
};

struct FlattenedSceneData
{
    std::vector<OptixInstance> instances;
    std::vector<LaunchParams::InstanceGeomRef> refs;
    std::vector<CUdeviceptr> ownedBuffers;
    ybi::float3 boundsMin = ybi::make_float3(std::numeric_limits<float>::max());
    ybi::float3 boundsMax = ybi::make_float3(-std::numeric_limits<float>::max());
};

static constexpr unsigned int kHitgroupSbtOffset = 0u;

static ybi::float3x4 IdentityTransform3x4()
{
    return ybi::float3x4(1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f);
}

static ybi::float4x4 ToFloat4x4(const ybi::float3x4 &m)
{
    return ybi::float4x4(m.m[0][0],
                         m.m[0][1],
                         m.m[0][2],
                         m.m[0][3],
                         m.m[1][0],
                         m.m[1][1],
                         m.m[1][2],
                         m.m[1][3],
                         m.m[2][0],
                         m.m[2][1],
                         m.m[2][2],
                         m.m[2][3],
                         0.0f,
                         0.0f,
                         0.0f,
                         1.0f);
}

static ybi::float3x4 ToFloat3x4(const ybi::float4x4 &m)
{
    return ybi::float3x4(m.m[0][0],
                         m.m[0][1],
                         m.m[0][2],
                         m.m[0][3],
                         m.m[1][0],
                         m.m[1][1],
                         m.m[1][2],
                         m.m[1][3],
                         m.m[2][0],
                         m.m[2][1],
                         m.m[2][2],
                         m.m[2][3]);
}

static void CopyTransform(const ybi::float3x4 &transform, float dst[12])
{
    std::memcpy(dst, transform.m, sizeof(float) * 12);
}

static void SetInstanceDefaults(OptixInstance &instance)
{
    std::memset(&instance, 0, sizeof(instance));
    instance.visibilityMask = 0xFFu;
    instance.flags = OPTIX_INSTANCE_FLAG_NONE;
}

static ybi::float3 TransformPoint(const ybi::float3x4 &m, const ybi::float3 &p)
{
    const ybi::float4 hp = ybi::make_float4(p.x, p.y, p.z, 1.0f);
    return ybi::mul(m, hp);
}

static void ExpandBounds(ybi::float3 &boundsMin, ybi::float3 &boundsMax, const ybi::float3 &p)
{
    boundsMin.x = std::min(boundsMin.x, p.x);
    boundsMin.y = std::min(boundsMin.y, p.y);
    boundsMin.z = std::min(boundsMin.z, p.z);
    boundsMax.x = std::max(boundsMax.x, p.x);
    boundsMax.y = std::max(boundsMax.y, p.y);
    boundsMax.z = std::max(boundsMax.z, p.z);
}

static void ExpandTransformedBounds(ybi::float3 &boundsMin,
                                    ybi::float3 &boundsMax,
                                    const ybi::float3x4 &worldFromLocal,
                                    const ybi::float3 &localMin,
                                    const ybi::float3 &localMax)
{
    const ybi::float3 corners[8] = {
        ybi::make_float3(localMin.x, localMin.y, localMin.z),
        ybi::make_float3(localMin.x, localMin.y, localMax.z),
        ybi::make_float3(localMin.x, localMax.y, localMin.z),
        ybi::make_float3(localMin.x, localMax.y, localMax.z),
        ybi::make_float3(localMax.x, localMin.y, localMin.z),
        ybi::make_float3(localMax.x, localMin.y, localMax.z),
        ybi::make_float3(localMax.x, localMax.y, localMin.z),
        ybi::make_float3(localMax.x, localMax.y, localMax.z),
    };
    for (const ybi::float3 &corner : corners)
    {
        ExpandBounds(boundsMin, boundsMax, TransformPoint(worldFromLocal, corner));
    }
}

static OptixTraversableHandle BuildTopLevelIAS(CUDADevice *device,
                                               HostMemoryArena &hostArena,
                                               const std::vector<OptixInstance> &instances)
{
    YBI_ASSERT(!instances.empty());
    CUDA_ASSERT(cuCtxPushCurrent(device->cudaContext));

    DeviceMemoryView<OptixInstance> deviceInstances =
        device->deviceArena->PushArray<OptixInstance>(instances.size());
    CUDA_ASSERT(cuMemcpyHtoD(
        CUdeviceptr(deviceInstances.data()), instances.data(), deviceInstances.numBytes()));

    OptixBuildInput buildInput = {};
    buildInput.type = OPTIX_BUILD_INPUT_TYPE_INSTANCES;
    buildInput.instanceArray.instances = CUdeviceptr(deviceInstances.data());
    buildInput.instanceArray.numInstances = (unsigned int)instances.size();

    OptixAccelBuildOptions options = {};
    options.buildFlags = OPTIX_BUILD_FLAG_PREFER_FAST_TRACE;
    options.operation = OPTIX_BUILD_OPERATION_BUILD;
    options.motionOptions.numKeys = 1;
    options.motionOptions.flags = 0;
    options.motionOptions.timeBegin = 0.0f;
    options.motionOptions.timeEnd = 1.0f;

    OptixAccelBufferSizes sizes = {};
    OPTIX_CHECK(optixAccelComputeMemoryUsage(
        device->optixDeviceContext, &options, &buildInput, 1, &sizes));

    DeviceMemoryView<uint8_t> tempBuffer = device->deviceArena->PushArray<uint8_t>(
        sizes.tempSizeInBytes, OPTIX_ACCEL_BUFFER_BYTE_ALIGNMENT);
    DeviceMemoryView<uint8_t> outputBuffer = device->Alloc<uint8_t>(sizes.outputSizeInBytes);

    OptixTraversableHandle outputHandle = {};
    OPTIX_CHECK(optixAccelBuild(device->optixDeviceContext,
                                0,
                                &options,
                                &buildInput,
                                1,
                                CUdeviceptr(tempBuffer.data()),
                                sizes.tempSizeInBytes,
                                CUdeviceptr(outputBuffer.data()),
                                sizes.outputSizeInBytes,
                                &outputHandle,
                                nullptr,
                                0));
    CUDA_ASSERT(cuStreamSynchronize(0));
    CUDA_ASSERT(cuCtxPopCurrent(0));
    hostArena.Clear();
    device->deviceArena->Clear();
    return outputHandle;
}

static void
AppendFlattenedSceneInstances(Scene *scene,
                              const ybi::float4x4 &worldFromScene,
                              const std::unordered_map<Scene *, size_t> &sceneIndexMap,
                              const std::vector<std::vector<MeshAccelData>> &sceneMeshAccels,
                              const std::vector<std::vector<CurveAccelData>> &sceneCurveAccels,
                              FlattenedSceneData &out)
{
    const auto sceneIndexIt = sceneIndexMap.find(scene);
    YBI_ASSERT(sceneIndexIt != sceneIndexMap.end());
    const size_t sceneIndex = sceneIndexIt->second;
    const std::vector<MeshAccelData> &meshAccels = sceneMeshAccels[sceneIndex];
    const std::vector<CurveAccelData> &curveAccels = sceneCurveAccels[sceneIndex];
    YBI_ASSERT(meshAccels.size() == scene->meshes.size());
    YBI_ASSERT(curveAccels.size() == scene->curves.size());

    for (size_t meshIndex = 0; meshIndex < scene->meshes.size(); meshIndex++)
    {
        const MeshAccelData &meshAccel = meshAccels[meshIndex];
        if (!meshAccel.gasHandle)
        {
            continue;
        }

        const ybi::float4x4 worldFromMesh =
            ybi::mul(worldFromScene, ToFloat4x4(scene->meshes[meshIndex].parentFromLocal));
        const ybi::float3x4 worldFromMesh3x4 = ToFloat3x4(worldFromMesh);

        OptixInstance optixInstance = {};
        SetInstanceDefaults(optixInstance);
        CopyTransform(worldFromMesh3x4, optixInstance.transform);
        optixInstance.instanceId = (unsigned int)out.refs.size();
        YBI_ASSERT(optixInstance.instanceId < (1u << 24));
        optixInstance.sbtOffset = kHitgroupSbtOffset;
        optixInstance.traversableHandle = meshAccel.gasHandle;
        out.instances.push_back(optixInstance);
        out.refs.push_back({(unsigned long long)meshAccel.positionsBuffer,
                            (unsigned long long)meshAccel.indicesBuffer,
                            meshAccel.numPositions,
                            meshAccel.numIndices});

        ExpandTransformedBounds(out.boundsMin,
                                out.boundsMax,
                                worldFromMesh3x4,
                                meshAccel.boundsMin,
                                meshAccel.boundsMax);
    }

    for (size_t curveIndex = 0; curveIndex < scene->curves.size(); curveIndex++)
    {
        const CurveAccelData &curveAccel = curveAccels[curveIndex];
        if (!curveAccel.gasHandle)
        {
            continue;
        }

        const ybi::float4x4 worldFromCurve =
            ybi::mul(worldFromScene, ToFloat4x4(scene->curves[curveIndex].parentFromLocal));
        const ybi::float3x4 worldFromCurve3x4 = ToFloat3x4(worldFromCurve);

        OptixInstance optixInstance = {};
        SetInstanceDefaults(optixInstance);
        CopyTransform(worldFromCurve3x4, optixInstance.transform);
        optixInstance.instanceId = (unsigned int)out.refs.size();
        YBI_ASSERT(optixInstance.instanceId < (1u << 24));
        optixInstance.sbtOffset = kHitgroupSbtOffset;
        optixInstance.traversableHandle = curveAccel.gasHandle;
        out.instances.push_back(optixInstance);
        out.refs.push_back({0ull, 0ull, 0, 0});

        ExpandTransformedBounds(out.boundsMin,
                                out.boundsMax,
                                worldFromCurve3x4,
                                curveAccel.boundsMin,
                                curveAccel.boundsMax);
    }

    for (const Instance &instance : scene->instances)
    {
        YBI_ASSERT(instance.childSceneIndex < scene->childScenes.size());
        Scene *childScene = scene->childScenes[instance.childSceneIndex];
        YBI_ASSERT(childScene != nullptr);
        const ybi::float4x4 worldFromChild =
            ybi::mul(worldFromScene, ToFloat4x4(instance.parentFromLocal));
        AppendFlattenedSceneInstances(
            childScene, worldFromChild, sceneIndexMap, sceneMeshAccels, sceneCurveAccels, out);
    }
}

static FlattenedSceneData
BuildFlattenedUSDScene(CUDADevice *device, HostMemoryArena &hostArena, ScenePool *scenePool)
{
    YBI_ASSERT(scenePool);
    YBI_ASSERT(scenePool->rootSceneIndex < scenePool->scenes.size());

    std::unordered_map<Scene *, size_t> sceneIndexMap;
    sceneIndexMap.reserve(scenePool->scenes.size());
    std::vector<std::vector<MeshAccelData>> sceneMeshAccels(scenePool->scenes.size());
    std::vector<std::vector<CurveAccelData>> sceneCurveAccels(scenePool->scenes.size());

    for (size_t sceneIndex = 0; sceneIndex < scenePool->scenes.size(); sceneIndex++)
    {
        Scene *scene = scenePool->scenes[sceneIndex].get();
        YBI_ASSERT(scene != nullptr);
        sceneIndexMap[scene] = sceneIndex;

        std::vector<MeshAccelData> &meshAccels = sceneMeshAccels[sceneIndex];
        meshAccels.resize(scene->meshes.size());
        for (size_t meshIndex = 0; meshIndex < scene->meshes.size(); meshIndex++)
        {
            Mesh &mesh = scene->meshes[meshIndex];
            if (mesh.positions.size() == 0 || mesh.indices.size() == 0)
            {
                continue;
            }

            std::vector<ybi::float3> meshPositions(mesh.positions.begin(), mesh.positions.end());
            std::vector<int> meshIndices(mesh.indices.begin(), mesh.indices.end());
            Mesh localMesh{Array<ybi::float3>(meshPositions), Array<int>(meshIndices)};
            MeshAccelData &meshAccel = meshAccels[meshIndex];
            meshAccel.gasHandle = BuildTriangleGASFromMesh(device, hostArena, localMesh);
            meshAccel.numPositions = (int)mesh.positions.size();
            meshAccel.numIndices = (int)mesh.indices.size();
            ComputeBounds(localMesh, meshAccel.boundsMin, meshAccel.boundsMax);

            CUDA_ASSERT(cuMemAlloc(&meshAccel.positionsBuffer,
                                   sizeof(ybi::float3) * mesh.positions.size()));
            CUDA_ASSERT(cuMemAlloc(&meshAccel.indicesBuffer, sizeof(int) * mesh.indices.size()));
            CUDA_ASSERT(cuMemcpyHtoD(meshAccel.positionsBuffer,
                                     mesh.positions.data(),
                                     sizeof(ybi::float3) * mesh.positions.size()));
            CUDA_ASSERT(cuMemcpyHtoD(
                meshAccel.indicesBuffer, mesh.indices.data(), sizeof(int) * mesh.indices.size()));

            hostArena.Clear();
            device->deviceArena->Clear();
        }

        std::vector<CurveAccelData> &curveAccels = sceneCurveAccels[sceneIndex];
        curveAccels.resize(scene->curves.size());
        for (size_t curveIndex = 0; curveIndex < scene->curves.size(); curveIndex++)
        {
            Curves &curves = scene->curves[curveIndex];
            if (curves.GetNumVertices() == 0 || curves.GetNumCurves() == 0)
            {
                continue;
            }

            const Array<ybi::float3> &curvePositions = curves.GetVertices();
            const Array<float> &curveWidths = curves.GetWidths();
            std::vector<ybi::float3> localPositions(curvePositions.begin(), curvePositions.end());
            std::vector<float> localWidths(curveWidths.begin(), curveWidths.end());
            std::vector<int> offsets;
            offsets.reserve(curves.GetNumCurves());
            for (size_t i = 0; i < curves.GetNumCurves(); i++)
            {
                offsets.push_back(curves.GetCurveKeyStart(i));
            }

            Curves localCurves{Array<ybi::float3>(localPositions),
                               Array<float>(localWidths),
                               Array<int>(offsets)};

            CurveAccelData &curveAccel = curveAccels[curveIndex];
            curveAccel.gasHandle = BuildCurveGASFromCurves(device, hostArena, localCurves);
            ComputeBounds(curves, curveAccel.boundsMin, curveAccel.boundsMax);
            hostArena.Clear();
            device->deviceArena->Clear();
        }
    }

    FlattenedSceneData result = {};
    Scene *rootScene = scenePool->scenes[scenePool->rootSceneIndex].get();
    AppendFlattenedSceneInstances(rootScene,
                                  ybi::float4x4::Identity(),
                                  sceneIndexMap,
                                  sceneMeshAccels,
                                  sceneCurveAccels,
                                  result);
    for (const std::vector<MeshAccelData> &meshAccels : sceneMeshAccels)
    {
        for (const MeshAccelData &meshAccel : meshAccels)
        {
            if (meshAccel.positionsBuffer)
            {
                result.ownedBuffers.push_back(meshAccel.positionsBuffer);
            }
            if (meshAccel.indicesBuffer)
            {
                result.ownedBuffers.push_back(meshAccel.indicesBuffer);
            }
        }
    }
    return result;
}

static void RenderTraversable(OptixPipeline pipeline,
                              const OptixShaderBindingTable &sbt,
                              OptixTraversableHandle traversable,
                              const ybi::float3 &boundsMin,
                              const ybi::float3 &boundsMax,
                              const char *outputFile,
                              IntegratorType integrator,
                              int spp,
                              CUdeviceptr instanceGeomRefsBuffer,
                              int instanceGeomRefCount,
                              const std::optional<RenderCameraOverride> &cameraOverride,
                              const std::optional<ybi::float3> &cameraPositionOverride,
                              const std::optional<ybi::float3> &lookAtOverride)
{
    printf("render: begin\n");
    fflush(stdout);
    const int width = cameraOverride.has_value() ? cameraOverride->width : 1280;
    const int height = cameraOverride.has_value() ? cameraOverride->height : 720;
    const size_t imageSize = static_cast<size_t>(width) * static_cast<size_t>(height) * 4;
    CUdeviceptr imageBuffer = 0;
    CUDA_ASSERT(cuMemAlloc(&imageBuffer, imageSize));
    printf("render: image buffer allocated\n");
    fflush(stdout);

    const ybi::float3 center = (boundsMin + boundsMax) * 0.5f;
    const ybi::float3 extent = boundsMax - boundsMin;
    const float diagonal = std::max(0.001f, ybi::length(extent));

    LaunchParams params = {};
    params.traversable = traversable;
    params.image = imageBuffer;
    params.width = width;
    params.height = height;
    if (cameraOverride.has_value())
    {
        params.cameraOrigin = cameraOverride->origin;
        params.cameraU = cameraOverride->U;
        params.cameraV = cameraOverride->V;
        params.cameraW = cameraOverride->W;
    }
    else
    {
        const ybi::float3 eye = cameraPositionOverride.has_value()
                                    ? cameraPositionOverride.value()
                                    : center + ybi::make_float3(0.0f, 0.0f, 1.25f * diagonal);
        const ybi::float3 lookAt = lookAtOverride.has_value() ? lookAtOverride.value() : center;
        const ybi::float3 forward = Normalize(lookAt - eye);
        ybi::float3 worldUp = ybi::make_float3(0.0f, 0.0f, 1.0f);
        if (std::abs(ybi::dot(forward, worldUp)) > 0.999f)
        {
            worldUp = ybi::make_float3(0.0f, 1.0f, 0.0f);
        }
        const ybi::float3 right = Normalize(Cross(forward, worldUp));
        const ybi::float3 up = Normalize(Cross(right, forward));
        const float aspect = static_cast<float>(width) / static_cast<float>(height);
        const float fovY = 45.0f * 3.14159265358979323846f / 180.0f;
        const float tanHalfFov = std::tan(fovY * 0.5f);
        params.cameraOrigin = eye;
        params.cameraU = right * (aspect * tanHalfFov);
        params.cameraV = up * tanHalfFov;
        params.cameraW = forward;
    }
    params.wireframe.lineWidth = 0.012f;
    params.wireframe.lineFeather = 0.006f;
    params.wireframe.edgeDarkness = 0.10f;
    params.wireframe.padding = 0.0f;
    params.integrator = integrator == IntegratorType::AO ? 1 : 0;
    params.spp = std::max(1, spp);
    params.aoBias = 0.002f * diagonal;
    params.aoMaxDistance = 0.25f * diagonal;
    params.instanceGeomRefs = (unsigned long long)instanceGeomRefsBuffer;
    params.instanceGeomRefCount = instanceGeomRefCount;

    CUdeviceptr paramsBuffer = 0;
    CUDA_ASSERT(cuMemAlloc(&paramsBuffer, sizeof(LaunchParams)));
    CUDA_ASSERT(cuMemcpyHtoD(paramsBuffer, &params, sizeof(LaunchParams)));
    printf("render: params uploaded\n");
    fflush(stdout);

    OPTIX_CHECK(
        optixLaunch(pipeline, 0, paramsBuffer, sizeof(LaunchParams), &sbt, width, height, 1));
    printf("render: launch returned\n");
    fflush(stdout);
    CUDA_ASSERT(cuStreamSynchronize(0));
    printf("render: stream synced\n");
    fflush(stdout);

    std::vector<uint8_t> hostImage(imageSize);
    CUDA_ASSERT(cuMemcpyDtoH(hostImage.data(), imageBuffer, imageSize));
    printf("render: copied image to host\n");
    fflush(stdout);
    const bool writeOk = SavePNG(outputFile, hostImage, width, height);
    printf("render: save returned=%d file=%s\n", writeOk ? 1 : 0, outputFile);
    fflush(stdout);
    if (!writeOk)
    {
        fprintf(stderr, "Failed to write PNG: %s\n", outputFile);
    }

    CUDA_ASSERT(cuMemFree(paramsBuffer));
    CUDA_ASSERT(cuMemFree(imageBuffer));
    printf("render: end\n");
    fflush(stdout);
}
} // namespace

int main(int argc, char **argv)
{
    printf("optix_harness: start\n");
    fflush(stdout);
    const CliOptions options = ParseCli(argc, argv);
    printf("optix_harness: parsed cli\n");
    fflush(stdout);

    CUDADevice device;
    printf("optix_harness: cuda device ready\n");
    fflush(stdout);
    HostMemoryArena hostArena;
    const std::string ptx = ReadTextFile(YBI_OPTIX_PRIMARY_PTX_PATH);
    printf("optix_harness: ptx loaded\n");
    fflush(stdout);

    if (!device.CreateOptixPrimaryPipeline(ptx))
    {
        fprintf(stderr, "Failed to create OptiX primary pipeline.\n");
        return 1;
    }
    const OptixPipeline pipeline = device.optixPrimaryPipeline.pipeline;
    const OptixShaderBindingTable &sbt = device.optixPrimaryPipeline.sbt;
    printf("optix_harness: pipeline created\n");
    fflush(stdout);

    if (options.type == RenderType::Triangle || options.type == RenderType::Cluster)
    {
        printf("optix_harness: loading obj %s\n", options.inputPath.c_str());
        fflush(stdout);
        Mesh mesh = LoadObjMesh(options.inputPath);
        printf("optix_harness: obj loaded verts=%zu tris=%zu\n",
               mesh.positions.size(),
               mesh.indices.size() / 3);
        fflush(stdout);
        if (mesh.positions.size() == 0 || mesh.indices.size() == 0)
        {
            fprintf(stderr, "OBJ has no geometry: %s\n", options.inputPath.c_str());
            return 1;
        }
        ybi::float3 meshBoundsMin, meshBoundsMax;
        ComputeBounds(mesh, meshBoundsMin, meshBoundsMax);

        CUdeviceptr meshPositionsBuffer = 0;
        CUdeviceptr meshIndicesBuffer = 0;
        CUdeviceptr instanceGeomRefsBuffer = 0;
        CUDA_ASSERT(cuMemAlloc(&meshPositionsBuffer, sizeof(ybi::float3) * mesh.positions.size()));
        CUDA_ASSERT(cuMemAlloc(&meshIndicesBuffer, sizeof(int) * mesh.indices.size()));
        CUDA_ASSERT(cuMemcpyHtoD(meshPositionsBuffer,
                                 mesh.positions.data(),
                                 sizeof(ybi::float3) * mesh.positions.size()));
        CUDA_ASSERT(cuMemcpyHtoD(
            meshIndicesBuffer, mesh.indices.data(), sizeof(int) * mesh.indices.size()));
        const LaunchParams::InstanceGeomRef meshRef = {
            (unsigned long long)meshPositionsBuffer,
            (unsigned long long)meshIndicesBuffer,
            (int)mesh.positions.size(),
            (int)mesh.indices.size(),
        };
        CUDA_ASSERT(cuMemAlloc(&instanceGeomRefsBuffer, sizeof(LaunchParams::InstanceGeomRef)));
        CUDA_ASSERT(
            cuMemcpyHtoD(instanceGeomRefsBuffer, &meshRef, sizeof(LaunchParams::InstanceGeomRef)));

        if (options.type == RenderType::Triangle)
        {
            printf("optix_harness: building triangle gas\n");
            fflush(stdout);
            OptixTraversableHandle triangleHandle =
                BuildTriangleGASFromMesh(&device, hostArena, mesh);
            if (triangleHandle)
            {
                OptixInstance rootInstance = {};
                SetInstanceDefaults(rootInstance);
                const ybi::float3x4 identity = IdentityTransform3x4();
                CopyTransform(identity, rootInstance.transform);
                rootInstance.instanceId = 0;
                rootInstance.sbtOffset = kHitgroupSbtOffset;
                rootInstance.traversableHandle = triangleHandle;

                const std::vector<OptixInstance> rootInstances = {rootInstance};
                OptixTraversableHandle rootIAS =
                    BuildTopLevelIAS(&device, hostArena, rootInstances);

                printf("optix_harness: rendering triangle ias\n");
                fflush(stdout);
                RenderTraversable(pipeline,
                                  sbt,
                                  rootIAS,
                                  meshBoundsMin,
                                  meshBoundsMax,
                                  options.outputPath.c_str(),
                                  options.integrator,
                                  options.spp,
                                  instanceGeomRefsBuffer,
                                  1,
                                  std::nullopt,
                                  options.cameraPosition,
                                  options.lookAt);
                printf("Wrote %s\n", options.outputPath.c_str());
            }
            else
            {
                printf("Triangle handle invalid; skipped triangle render.\n");
            }
        }
        else
        {
#if (OPTIX_VERSION >= 90000)
            printf("optix_harness: building cluster gas\n");
            fflush(stdout);
            OptixTraversableHandle clusterHandle =
                BuildClusterGASFromMesh(&device, hostArena, mesh);
            if (clusterHandle)
            {
                OptixInstance rootInstance = {};
                SetInstanceDefaults(rootInstance);
                const ybi::float3x4 identity = IdentityTransform3x4();
                CopyTransform(identity, rootInstance.transform);
                rootInstance.instanceId = 0;
                rootInstance.sbtOffset = kHitgroupSbtOffset;
                rootInstance.traversableHandle = clusterHandle;

                const std::vector<OptixInstance> rootInstances = {rootInstance};
                OptixTraversableHandle rootIAS =
                    BuildTopLevelIAS(&device, hostArena, rootInstances);

                printf("optix_harness: rendering cluster ias\n");
                fflush(stdout);
                RenderTraversable(pipeline,
                                  sbt,
                                  rootIAS,
                                  meshBoundsMin,
                                  meshBoundsMax,
                                  options.outputPath.c_str(),
                                  options.integrator,
                                  options.spp,
                                  instanceGeomRefsBuffer,
                                  1,
                                  std::nullopt,
                                  options.cameraPosition,
                                  options.lookAt);
                printf("Wrote %s\n", options.outputPath.c_str());
            }
            else
            {
                printf("Cluster handle invalid; skipped cluster render.\n");
            }
#else
            printf("Cluster GAS not supported on this OptiX version.\n");
#endif
        }

        CUDA_ASSERT(cuMemFree(instanceGeomRefsBuffer));
        CUDA_ASSERT(cuMemFree(meshIndicesBuffer));
        CUDA_ASSERT(cuMemFree(meshPositionsBuffer));
    }
    else if (options.type == RenderType::Curve)
    {
        printf("optix_harness: loading curve json %s\n", options.inputPath.c_str());
        fflush(stdout);
        Curves curves = LoadCurveJson(options.inputPath);
        printf("optix_harness: curve loaded verts=%zu curves=%zu\n",
               curves.GetNumVertices(),
               curves.GetNumCurves());
        fflush(stdout);
        if (curves.GetNumVertices() > 0 && curves.GetNumCurves() > 0)
        {
            ybi::float3 curveBoundsMin, curveBoundsMax;
            ComputeBounds(curves, curveBoundsMin, curveBoundsMax);
            OptixTraversableHandle curveHandle =
                BuildCurveGASFromCurves(&device, hostArena, curves);
            if (curveHandle)
            {
                OptixInstance rootInstance = {};
                SetInstanceDefaults(rootInstance);
                const ybi::float3x4 identity = IdentityTransform3x4();
                CopyTransform(identity, rootInstance.transform);
                rootInstance.instanceId = 0;
                rootInstance.sbtOffset = kHitgroupSbtOffset;
                rootInstance.traversableHandle = curveHandle;

                const std::vector<OptixInstance> rootInstances = {rootInstance};
                OptixTraversableHandle rootIAS =
                    BuildTopLevelIAS(&device, hostArena, rootInstances);

                RenderTraversable(pipeline,
                                  sbt,
                                  rootIAS,
                                  curveBoundsMin,
                                  curveBoundsMax,
                                  options.outputPath.c_str(),
                                  options.integrator,
                                  options.spp,
                                  0,
                                  0,
                                  std::nullopt,
                                  options.cameraPosition,
                                  options.lookAt);
                printf("Wrote %s\n", options.outputPath.c_str());
            }
            else
            {
                printf("Curve handle invalid; skipped curve render.\n");
            }
        }
        else
        {
            printf("Curve JSON empty or invalid; skipped curve render: %s\n",
                   options.inputPath.c_str());
        }
    }
    else
    {
        printf("optix_harness: loading usd scene %s\n", options.inputPath.c_str());
        fflush(stdout);

        ScenePool scenePool = {};
        LoadUSDScene(&scenePool, options.inputPath);
        if (scenePool.scenes.empty() || scenePool.rootSceneIndex >= scenePool.scenes.size())
        {
            fprintf(stderr,
                    "Failed to load USD scene or invalid root: %s\n",
                    options.inputPath.c_str());
            return 1;
        }

        FlattenedSceneData flattened = BuildFlattenedUSDScene(&device, hostArena, &scenePool);
        if (flattened.instances.empty() || flattened.refs.empty())
        {
            fprintf(
                stderr, "USD scene produced no mesh instances: %s\n", options.inputPath.c_str());
            return 1;
        }
        YBI_ASSERT(flattened.instances.size() == flattened.refs.size());

        OptixTraversableHandle rootIAS = BuildTopLevelIAS(&device, hostArena, flattened.instances);

        CUdeviceptr instanceGeomRefsBuffer = 0;
        CUDA_ASSERT(cuMemAlloc(&instanceGeomRefsBuffer,
                               flattened.refs.size() * sizeof(LaunchParams::InstanceGeomRef)));
        CUDA_ASSERT(cuMemcpyHtoD(instanceGeomRefsBuffer,
                                 flattened.refs.data(),
                                 flattened.refs.size() * sizeof(LaunchParams::InstanceGeomRef)));

        std::optional<RenderCameraOverride> usdCamera = std::nullopt;
        if (!options.cameraPosition.has_value() && !options.lookAt.has_value())
        {
            usdCamera = BuildUsdRenderCamera(scenePool.camera);
            if (usdCamera.has_value())
            {
                printf("optix_harness: using usd camera viewport=%dx%d\n",
                       usdCamera->width,
                       usdCamera->height);
            }
            else
            {
                printf("optix_harness: usd camera unavailable, using bounds camera\n");
            }
            fflush(stdout);
        }

        RenderTraversable(pipeline,
                          sbt,
                          rootIAS,
                          flattened.boundsMin,
                          flattened.boundsMax,
                          options.outputPath.c_str(),
                          options.integrator,
                          options.spp,
                          instanceGeomRefsBuffer,
                          (int)flattened.refs.size(),
                          usdCamera,
                          options.cameraPosition,
                          options.lookAt);
        printf("Wrote %s\n", options.outputPath.c_str());
        CUDA_ASSERT(cuMemFree(instanceGeomRefsBuffer));
        for (CUdeviceptr buffer : flattened.ownedBuffers)
        {
            CUDA_ASSERT(cuMemFree(buffer));
        }
    }
    hostArena.Clear();
    device.deviceArena->Clear();

    device.DestroyOptixPrimaryPipeline();
    std::fflush(stdout);
    std::fflush(stderr);
    std::_Exit(0);
}
