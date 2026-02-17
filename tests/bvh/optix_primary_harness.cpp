#include "device/cuda_device.h"
#include "scene/scene.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "third_party/stb_image_write.h"
#include "util/array.h"
#include "util/float3.h"
#include <algorithm>
#include <cctype>
#include <cstddef>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <limits>
#include <optional>
#include <sstream>
#include <string>
#include <vector>
#include <cstdlib>

#include <optix_stack_size.h>
#include <optix_stubs.h>

using namespace ybi;

#ifndef YBI_OPTIX_PRIMARY_PTX_PATH
#define YBI_OPTIX_PRIMARY_PTX_PATH ""
#endif

namespace
{
struct LaunchParams
{
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
};

enum class RenderType
{
    Triangle,
    Cluster,
    Curve
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

template <typename T>
struct alignas(OPTIX_SBT_RECORD_ALIGNMENT) SbtRecord
{
    char header[OPTIX_SBT_RECORD_HEADER_SIZE];
    T data;
};

struct EmptyData
{
};

struct HitgroupData
{
    unsigned long long positions;
    unsigned long long indices;
    int numPositions;
    int numIndices;
};

using RaygenRecord = SbtRecord<EmptyData>;
using MissRecord = SbtRecord<EmptyData>;
using HitgroupRecord = SbtRecord<HitgroupData>;

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
    return ybi::make_float3(
        a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
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
    printf("Usage: %s [--type triangle|cluster|curve] [--file path] [--out path] "
           "[--integrator primary|ao] [--spp N] [--cam-pos x y z] [--look-at x y z]\n",
           exeName);
    printf("  --type triangle|cluster|curve\n");
    printf("  --file OBJ path for triangle/cluster; JSON path for curve\n");
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

static void UploadHitgroupData(CUdeviceptr hitgroupRecordBuffer,
                               CUdeviceptr positionsBuffer,
                               CUdeviceptr indicesBuffer,
                               int numPositions,
                               int numIndices)
{
    HitgroupData hitgroupData = {};
    hitgroupData.positions = (unsigned long long)positionsBuffer;
    hitgroupData.indices = (unsigned long long)indicesBuffer;
    hitgroupData.numPositions = numPositions;
    hitgroupData.numIndices = numIndices;
    CUDA_ASSERT(cuMemcpyHtoD(hitgroupRecordBuffer + offsetof(HitgroupRecord, data),
                             &hitgroupData,
                             sizeof(HitgroupData)));
}

static bool ParseObjIndex(const std::string &token, int &indexOut)
{
    const size_t slashPos = token.find('/');
    const std::string indexText = slashPos == std::string::npos ? token : token.substr(0, slashPos);
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

    const std::string json(
        (std::istreambuf_iterator<char>(input)), std::istreambuf_iterator<char>());
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
        positions.push_back(ybi::make_float3(pointScalars[i + 0], pointScalars[i + 1], pointScalars[i + 2]));
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

static OptixPipeline CreatePipeline(OptixDeviceContext optixContext,
                                    const std::string &ptx,
                                    bool useCurveIntersection,
                                    OptixShaderBindingTable &sbtOut,
                                    CUdeviceptr &raygenRecordBufferOut,
                                    CUdeviceptr &missRecordBufferOut,
                                    CUdeviceptr &hitgroupRecordBufferOut,
                                    OptixModule &moduleOut,
                                    OptixModule &curveModuleOut,
                                    OptixProgramGroup &raygenGroupOut,
                                    OptixProgramGroup &missGroupOut,
                                    OptixProgramGroup &hitgroupGroupOut)
{
    curveModuleOut = nullptr;
    OptixModuleCompileOptions moduleCompileOptions = {};
    moduleCompileOptions.maxRegisterCount = OPTIX_COMPILE_DEFAULT_MAX_REGISTER_COUNT;
    moduleCompileOptions.optLevel = OPTIX_COMPILE_OPTIMIZATION_DEFAULT;
    moduleCompileOptions.debugLevel = OPTIX_COMPILE_DEBUG_LEVEL_DEFAULT;

    OptixPipelineCompileOptions pipelineCompileOptions = {};
    pipelineCompileOptions.traversableGraphFlags = OPTIX_TRAVERSABLE_GRAPH_FLAG_ALLOW_SINGLE_GAS;
    pipelineCompileOptions.usesMotionBlur = 0;
    pipelineCompileOptions.usesPrimitiveTypeFlags = useCurveIntersection
                                                        ? (OPTIX_PRIMITIVE_TYPE_FLAGS_ROUND_LINEAR |
                                                           OPTIX_PRIMITIVE_TYPE_FLAGS_ROUND_CUBIC_BSPLINE)
                                                        : OPTIX_PRIMITIVE_TYPE_FLAGS_TRIANGLE;
    pipelineCompileOptions.numPayloadValues = 1;
    pipelineCompileOptions.numAttributeValues = useCurveIntersection ? 4 : 2;
    pipelineCompileOptions.exceptionFlags = OPTIX_EXCEPTION_FLAG_NONE;
    pipelineCompileOptions.pipelineLaunchParamsVariableName = "params";

    char log[2048];
    size_t logSize = sizeof(log);
    OPTIX_ASSERT(optixModuleCreate(optixContext,
                                   &moduleCompileOptions,
                                   &pipelineCompileOptions,
                                   ptx.c_str(),
                                   ptx.size(),
                                   log,
                                   &logSize,
                                   &moduleOut));

    if (useCurveIntersection)
    {
        OptixBuiltinISOptions builtinISOptions = {};
        builtinISOptions.builtinISModuleType = OPTIX_PRIMITIVE_TYPE_ROUND_CUBIC_BSPLINE;
        builtinISOptions.usesMotionBlur = 0;
        builtinISOptions.buildFlags = OPTIX_BUILD_FLAG_NONE;
        builtinISOptions.curveEndcapFlags = OPTIX_CURVE_ENDCAP_DEFAULT;
        OPTIX_ASSERT(optixBuiltinISModuleGet(optixContext,
                                             &moduleCompileOptions,
                                             &pipelineCompileOptions,
                                             &builtinISOptions,
                                             &curveModuleOut));
    }

    OptixProgramGroupOptions programGroupOptions = {};
    OptixProgramGroupDesc raygenDesc = {};
    raygenDesc.kind = OPTIX_PROGRAM_GROUP_KIND_RAYGEN;
    raygenDesc.raygen.module = moduleOut;
    raygenDesc.raygen.entryFunctionName = "__raygen__primary";
    logSize = sizeof(log);
    OPTIX_ASSERT(optixProgramGroupCreate(
        optixContext, &raygenDesc, 1, &programGroupOptions, log, &logSize, &raygenGroupOut));

    OptixProgramGroupDesc missDesc = {};
    missDesc.kind = OPTIX_PROGRAM_GROUP_KIND_MISS;
    missDesc.miss.module = moduleOut;
    missDesc.miss.entryFunctionName = "__miss__primary";
    logSize = sizeof(log);
    OPTIX_ASSERT(optixProgramGroupCreate(
        optixContext, &missDesc, 1, &programGroupOptions, log, &logSize, &missGroupOut));

    OptixProgramGroupDesc hitgroupDesc = {};
    hitgroupDesc.kind = OPTIX_PROGRAM_GROUP_KIND_HITGROUP;
    hitgroupDesc.hitgroup.moduleCH = moduleOut;
    hitgroupDesc.hitgroup.entryFunctionNameCH = "__closesthit__primary";
    hitgroupDesc.hitgroup.moduleAH = moduleOut;
    hitgroupDesc.hitgroup.entryFunctionNameAH = "__anyhit__primary";
    hitgroupDesc.hitgroup.moduleIS = useCurveIntersection ? curveModuleOut : nullptr;
    hitgroupDesc.hitgroup.entryFunctionNameIS = nullptr;
    logSize = sizeof(log);
    OPTIX_ASSERT(optixProgramGroupCreate(
        optixContext, &hitgroupDesc, 1, &programGroupOptions, log, &logSize, &hitgroupGroupOut));

    OptixProgramGroup groups[] = {raygenGroupOut, missGroupOut, hitgroupGroupOut};
    OptixPipelineLinkOptions pipelineLinkOptions = {};
    pipelineLinkOptions.maxTraceDepth = 2;
    OptixPipeline pipeline = nullptr;
    logSize = sizeof(log);
    OPTIX_ASSERT(optixPipelineCreate(optixContext,
                                     &pipelineCompileOptions,
                                     &pipelineLinkOptions,
                                     groups,
                                     3,
                                     log,
                                     &logSize,
                                     &pipeline));

    OptixStackSizes stackSizes = {};
    OPTIX_ASSERT(optixUtilAccumulateStackSizes(raygenGroupOut, &stackSizes, pipeline));
    OPTIX_ASSERT(optixUtilAccumulateStackSizes(missGroupOut, &stackSizes, pipeline));
    OPTIX_ASSERT(optixUtilAccumulateStackSizes(hitgroupGroupOut, &stackSizes, pipeline));
    uint32_t directCallableStackSizeFromTraversal = 0;
    uint32_t directCallableStackSizeFromState = 0;
    uint32_t continuationStackSize = 0;
    OPTIX_ASSERT(optixUtilComputeStackSizes(&stackSizes,
                                            2,
                                            0,
                                            0,
                                            &directCallableStackSizeFromTraversal,
                                            &directCallableStackSizeFromState,
                                            &continuationStackSize));
    OPTIX_ASSERT(optixPipelineSetStackSize(pipeline,
                                           directCallableStackSizeFromTraversal,
                                           directCallableStackSizeFromState,
                                           continuationStackSize,
                                           1));

    RaygenRecord raygenRecord = {};
    MissRecord missRecord = {};
    HitgroupRecord hitgroupRecord = {};
    OPTIX_ASSERT(optixSbtRecordPackHeader(raygenGroupOut, &raygenRecord));
    OPTIX_ASSERT(optixSbtRecordPackHeader(missGroupOut, &missRecord));
    OPTIX_ASSERT(optixSbtRecordPackHeader(hitgroupGroupOut, &hitgroupRecord));

    CUDA_ASSERT(cuMemAlloc(&raygenRecordBufferOut, sizeof(RaygenRecord)));
    CUDA_ASSERT(cuMemAlloc(&missRecordBufferOut, sizeof(MissRecord)));
    CUDA_ASSERT(cuMemAlloc(&hitgroupRecordBufferOut, sizeof(HitgroupRecord)));
    CUDA_ASSERT(cuMemcpyHtoD(raygenRecordBufferOut, &raygenRecord, sizeof(RaygenRecord)));
    CUDA_ASSERT(cuMemcpyHtoD(missRecordBufferOut, &missRecord, sizeof(MissRecord)));
    CUDA_ASSERT(cuMemcpyHtoD(hitgroupRecordBufferOut, &hitgroupRecord, sizeof(HitgroupRecord)));

    sbtOut = {};
    sbtOut.raygenRecord = raygenRecordBufferOut;
    sbtOut.missRecordBase = missRecordBufferOut;
    sbtOut.missRecordStrideInBytes = sizeof(MissRecord);
    sbtOut.missRecordCount = 1;
    sbtOut.hitgroupRecordBase = hitgroupRecordBufferOut;
    sbtOut.hitgroupRecordStrideInBytes = sizeof(HitgroupRecord);
    sbtOut.hitgroupRecordCount = 1;
    return pipeline;
}

static void RenderTraversable(OptixPipeline pipeline,
                              const OptixShaderBindingTable &sbt,
                              OptixTraversableHandle traversable,
                              const ybi::float3 &boundsMin,
                              const ybi::float3 &boundsMax,
                              const char *outputFile,
                              IntegratorType integrator,
                              int spp,
                              const std::optional<ybi::float3> &cameraPositionOverride,
                              const std::optional<ybi::float3> &lookAtOverride)
{
    printf("render: begin\n");
    fflush(stdout);
    const int width = 1280;
    const int height = 720;
    const size_t imageSize = static_cast<size_t>(width) * static_cast<size_t>(height) * 4;
    CUdeviceptr imageBuffer = 0;
    CUDA_ASSERT(cuMemAlloc(&imageBuffer, imageSize));
    printf("render: image buffer allocated\n");
    fflush(stdout);

    const ybi::float3 center = (boundsMin + boundsMax) * 0.5f;
    const ybi::float3 extent = boundsMax - boundsMin;
    const float diagonal = std::max(0.001f, ybi::length(extent));

    const ybi::float3 eye =
        cameraPositionOverride.has_value() ? cameraPositionOverride.value()
                                           : center + ybi::make_float3(0.0f, 0.0f, 1.25f * diagonal);
    const ybi::float3 lookAt = lookAtOverride.has_value() ? lookAtOverride.value() : center;
    const ybi::float3 forward = Normalize(lookAt - eye);
    const ybi::float3 worldUp = ybi::make_float3(0.0f, 1.0f, 0.0f);
    const ybi::float3 right = Normalize(Cross(forward, worldUp));
    const ybi::float3 up = Normalize(Cross(right, forward));
    const float aspect = static_cast<float>(width) / static_cast<float>(height);
    const float fovY = 45.0f * 3.14159265358979323846f / 180.0f;
    const float tanHalfFov = std::tan(fovY * 0.5f);

    LaunchParams params = {};
    params.traversable = traversable;
    params.image = imageBuffer;
    params.width = width;
    params.height = height;
    params.cameraOrigin = eye;
    params.cameraU = right * (aspect * tanHalfFov);
    params.cameraV = up * tanHalfFov;
    params.cameraW = forward;
    params.wireframe.lineWidth = 0.012f;
    params.wireframe.lineFeather = 0.006f;
    params.wireframe.edgeDarkness = 0.10f;
    params.wireframe.padding = 0.0f;
    params.integrator = integrator == IntegratorType::AO ? 1 : 0;
    params.spp = std::max(1, spp);
    params.aoBias = 0.002f * diagonal;
    params.aoMaxDistance = 0.25f * diagonal;

    CUdeviceptr paramsBuffer = 0;
    CUDA_ASSERT(cuMemAlloc(&paramsBuffer, sizeof(LaunchParams)));
    CUDA_ASSERT(cuMemcpyHtoD(paramsBuffer, &params, sizeof(LaunchParams)));
    printf("render: params uploaded\n");
    fflush(stdout);

    OPTIX_ASSERT(optixLaunch(
        pipeline, 0, paramsBuffer, sizeof(LaunchParams), &sbt, width, height, 1));
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

    OptixShaderBindingTable sbt = {};
    CUdeviceptr raygenRecordBuffer = 0;
    CUdeviceptr missRecordBuffer = 0;
    CUdeviceptr hitgroupRecordBuffer = 0;
    OptixModule module = nullptr;
    OptixModule curveModule = nullptr;
    OptixProgramGroup raygenGroup = nullptr;
    OptixProgramGroup missGroup = nullptr;
    OptixProgramGroup hitgroupGroup = nullptr;
    OptixPipeline pipeline = nullptr;

    OptixShaderBindingTable curveSbt = {};
    CUdeviceptr curveRaygenRecordBuffer = 0;
    CUdeviceptr curveMissRecordBuffer = 0;
    CUdeviceptr curveHitgroupRecordBuffer = 0;
    OptixModule curvePipelineModule = nullptr;
    OptixModule curveIsModule = nullptr;
    OptixProgramGroup curveRaygenGroup = nullptr;
    OptixProgramGroup curveMissGroup = nullptr;
    OptixProgramGroup curveHitgroupGroup = nullptr;
    OptixPipeline curvePipeline = nullptr;

    if (options.type == RenderType::Triangle || options.type == RenderType::Cluster)
    {
        pipeline = CreatePipeline(device.optixDeviceContext,
                                  ptx,
                                  false,
                                  sbt,
                                  raygenRecordBuffer,
                                  missRecordBuffer,
                                  hitgroupRecordBuffer,
                                  module,
                                  curveModule,
                                  raygenGroup,
                                  missGroup,
                                  hitgroupGroup);
        printf("optix_harness: triangle pipeline created\n");
        fflush(stdout);
    }
    else
    {
        curvePipeline = CreatePipeline(device.optixDeviceContext,
                                       ptx,
                                       true,
                                       curveSbt,
                                       curveRaygenRecordBuffer,
                                       curveMissRecordBuffer,
                                       curveHitgroupRecordBuffer,
                                       curvePipelineModule,
                                       curveIsModule,
                                       curveRaygenGroup,
                                       curveMissGroup,
                                       curveHitgroupGroup);
        printf("optix_harness: curve pipeline created\n");
        fflush(stdout);
    }

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
        CUDA_ASSERT(cuMemAlloc(&meshPositionsBuffer, sizeof(ybi::float3) * mesh.positions.size()));
        CUDA_ASSERT(cuMemAlloc(&meshIndicesBuffer, sizeof(int) * mesh.indices.size()));
        CUDA_ASSERT(cuMemcpyHtoD(meshPositionsBuffer,
                                 mesh.positions.data(),
                                 sizeof(ybi::float3) * mesh.positions.size()));
        CUDA_ASSERT(cuMemcpyHtoD(meshIndicesBuffer, mesh.indices.data(), sizeof(int) * mesh.indices.size()));
        UploadHitgroupData(hitgroupRecordBuffer,
                           meshPositionsBuffer,
                           meshIndicesBuffer,
                           (int)mesh.positions.size(),
                           (int)mesh.indices.size());

        if (options.type == RenderType::Triangle)
        {
            printf("optix_harness: building triangle gas\n");
            fflush(stdout);
            OptixTraversableHandle triangleHandle = BuildTriangleGASFromMesh(&device, hostArena, mesh);
            printf("optix_harness: rendering triangle gas\n");
            fflush(stdout);
            RenderTraversable(pipeline,
                              sbt,
                              triangleHandle,
                              meshBoundsMin,
                              meshBoundsMax,
                              options.outputPath.c_str(),
                              options.integrator,
                              options.spp,
                              options.cameraPosition,
                              options.lookAt);
            printf("Wrote %s\n", options.outputPath.c_str());
        }
        else
        {
#if (OPTIX_VERSION >= 90000)
            printf("optix_harness: building cluster gas\n");
            fflush(stdout);
            OptixTraversableHandle clusterHandle = BuildClusterGASFromMesh(&device, hostArena, mesh);
            if (clusterHandle)
            {
                printf("optix_harness: rendering cluster gas\n");
                fflush(stdout);
                RenderTraversable(pipeline,
                                  sbt,
                                  clusterHandle,
                                  meshBoundsMin,
                                  meshBoundsMax,
                                  options.outputPath.c_str(),
                                  options.integrator,
                                  options.spp,
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

        CUDA_ASSERT(cuMemFree(meshIndicesBuffer));
        CUDA_ASSERT(cuMemFree(meshPositionsBuffer));
    }
    else
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
            OptixTraversableHandle curveHandle = BuildCurveGASFromCurves(&device, hostArena, curves);
            if (curveHandle)
            {
                RenderTraversable(curvePipeline,
                                  curveSbt,
                                  curveHandle,
                                  curveBoundsMin,
                                  curveBoundsMax,
                                  options.outputPath.c_str(),
                                  options.integrator,
                                  options.spp,
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
            printf("Curve JSON empty or invalid; skipped curve render: %s\n", options.inputPath.c_str());
        }
    }
    hostArena.Clear();
    device.deviceArena->Clear();

    if (pipeline)
    {
        CUDA_ASSERT(cuMemFree(hitgroupRecordBuffer));
        CUDA_ASSERT(cuMemFree(missRecordBuffer));
        CUDA_ASSERT(cuMemFree(raygenRecordBuffer));
        OPTIX_ASSERT(optixPipelineDestroy(pipeline));
        OPTIX_ASSERT(optixProgramGroupDestroy(hitgroupGroup));
        OPTIX_ASSERT(optixProgramGroupDestroy(missGroup));
        OPTIX_ASSERT(optixProgramGroupDestroy(raygenGroup));
        if (curveModule)
        {
            OPTIX_ASSERT(optixModuleDestroy(curveModule));
        }
        OPTIX_ASSERT(optixModuleDestroy(module));
    }

    if (curvePipeline)
    {
        CUDA_ASSERT(cuMemFree(curveHitgroupRecordBuffer));
        CUDA_ASSERT(cuMemFree(curveMissRecordBuffer));
        CUDA_ASSERT(cuMemFree(curveRaygenRecordBuffer));
        OPTIX_ASSERT(optixPipelineDestroy(curvePipeline));
        OPTIX_ASSERT(optixProgramGroupDestroy(curveHitgroupGroup));
        OPTIX_ASSERT(optixProgramGroupDestroy(curveMissGroup));
        OPTIX_ASSERT(optixProgramGroupDestroy(curveRaygenGroup));
        if (curveIsModule)
        {
            OPTIX_ASSERT(optixModuleDestroy(curveIsModule));
        }
        OPTIX_ASSERT(optixModuleDestroy(curvePipelineModule));
    }
    std::fflush(stdout);
    std::fflush(stderr);
    std::_Exit(0);
}
