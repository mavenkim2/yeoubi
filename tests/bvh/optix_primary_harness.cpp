#include "device/cuda_device.h"
#include "scene/scene.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "third_party/stb_image_write.h"
#include "util/array.h"
#include "util/float3.h"
#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

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

using RaygenRecord = SbtRecord<EmptyData>;
using MissRecord = SbtRecord<EmptyData>;
using HitgroupRecord = SbtRecord<EmptyData>;

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

static bool SavePNG(const char *filePath, const std::vector<uint8_t> &rgba, int width, int height)
{
    const int strideInBytes = width * 4;
    return stbi_write_png(filePath, width, height, 4, rgba.data(), strideInBytes) != 0;
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
    pipelineLinkOptions.maxTraceDepth = 1;
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
                                            1,
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
                              const char *outputFile)
{
    const int width = 1280;
    const int height = 720;
    const size_t imageSize = static_cast<size_t>(width) * static_cast<size_t>(height) * 4;
    CUdeviceptr imageBuffer = 0;
    CUDA_ASSERT(cuMemAlloc(&imageBuffer, imageSize));

    const ybi::float3 center = (boundsMin + boundsMax) * 0.5f;
    const ybi::float3 extent = boundsMax - boundsMin;
    const float diagonal = std::max(0.001f, ybi::length(extent));

    const ybi::float3 eye = center + ybi::make_float3(0.0f, 0.0f, 1.25f * diagonal);
    const ybi::float3 forward = Normalize(center - eye);
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

    CUdeviceptr paramsBuffer = 0;
    CUDA_ASSERT(cuMemAlloc(&paramsBuffer, sizeof(LaunchParams)));
    CUDA_ASSERT(cuMemcpyHtoD(paramsBuffer, &params, sizeof(LaunchParams)));

    OPTIX_ASSERT(optixLaunch(
        pipeline, 0, paramsBuffer, sizeof(LaunchParams), &sbt, width, height, 1));
    CUDA_ASSERT(cuStreamSynchronize(0));

    std::vector<uint8_t> hostImage(imageSize);
    CUDA_ASSERT(cuMemcpyDtoH(hostImage.data(), imageBuffer, imageSize));
    if (!SavePNG(outputFile, hostImage, width, height))
    {
        fprintf(stderr, "Failed to write PNG: %s\n", outputFile);
    }

    CUDA_ASSERT(cuMemFree(paramsBuffer));
    CUDA_ASSERT(cuMemFree(imageBuffer));
}
} // namespace

int main()
{
    const std::string objPath = "tests/bvh/out/stoat_body_selected.obj";
    const std::string curveJsonPath = "tests/bvh/out/selected_curve.json";
    Mesh mesh = LoadObjMesh(objPath);
    if (mesh.positions.size() == 0 || mesh.indices.size() == 0)
    {
        fprintf(stderr, "OBJ has no geometry: %s\n", objPath.c_str());
        return 1;
    }
    Curves curves = LoadCurveJson(curveJsonPath);

    CUDADevice device;
    HostMemoryArena hostArena;
    const std::string ptx = ReadTextFile(YBI_OPTIX_PRIMARY_PTX_PATH);

    OptixShaderBindingTable sbt = {};
    CUdeviceptr raygenRecordBuffer = 0;
    CUdeviceptr missRecordBuffer = 0;
    CUdeviceptr hitgroupRecordBuffer = 0;
    OptixModule module = nullptr;
    OptixModule curveModule = nullptr;
    OptixProgramGroup raygenGroup = nullptr;
    OptixProgramGroup missGroup = nullptr;
    OptixProgramGroup hitgroupGroup = nullptr;
    OptixPipeline pipeline = CreatePipeline(device.optixDeviceContext,
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

    OptixShaderBindingTable curveSbt = {};
    CUdeviceptr curveRaygenRecordBuffer = 0;
    CUdeviceptr curveMissRecordBuffer = 0;
    CUdeviceptr curveHitgroupRecordBuffer = 0;
    OptixModule curvePipelineModule = nullptr;
    OptixModule curveIsModule = nullptr;
    OptixProgramGroup curveRaygenGroup = nullptr;
    OptixProgramGroup curveMissGroup = nullptr;
    OptixProgramGroup curveHitgroupGroup = nullptr;
    OptixPipeline curvePipeline = CreatePipeline(device.optixDeviceContext,
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

    ybi::float3 meshBoundsMin, meshBoundsMax;
    ComputeBounds(mesh, meshBoundsMin, meshBoundsMax);

    OptixTraversableHandle triangleHandle = BuildTriangleGASFromMesh(&device, hostArena, mesh);
    RenderTraversable(
        pipeline, sbt, triangleHandle, meshBoundsMin, meshBoundsMax, "optix_triangle_gas.png");
    printf("Wrote optix_triangle_gas.png\n");
    hostArena.Clear();
    device.deviceArena->Clear();

#if (OPTIX_VERSION >= 90000)
    OptixTraversableHandle clusterHandle = BuildClusterGASFromMesh(&device, hostArena, mesh);
    if (clusterHandle)
    {
        RenderTraversable(
            pipeline, sbt, clusterHandle, meshBoundsMin, meshBoundsMax, "optix_cluster_gas.png");
        printf("Wrote optix_cluster_gas.png\n");
    }
    else
    {
        printf("Cluster handle invalid; skipped cluster render.\n");
    }
    hostArena.Clear();
    device.deviceArena->Clear();
#else
    printf("Cluster GAS not supported on this OptiX version.\n");
#endif

    if (curves.GetNumVertices() > 0 && curves.GetNumCurves() > 0)
    {
        ybi::float3 curveBoundsMin, curveBoundsMax;
        ComputeBounds(curves, curveBoundsMin, curveBoundsMax);
        printf("Curve bounds min=(%f,%f,%f) max=(%f,%f,%f)\n",
               curveBoundsMin.x,
               curveBoundsMin.y,
               curveBoundsMin.z,
               curveBoundsMax.x,
               curveBoundsMax.y,
               curveBoundsMax.z);
        OptixTraversableHandle curveHandle = BuildCurveGASFromCurves(&device, hostArena, curves);
        printf("Curve handle=%llu\n", (unsigned long long)curveHandle);
        if (curveHandle)
        {
            RenderTraversable(curvePipeline,
                              curveSbt,
                              curveHandle,
                              curveBoundsMin,
                              curveBoundsMax,
                              "optix_curve_gas.png");
            printf("Wrote optix_curve_gas.png\n");
        }
        else
        {
            printf("Curve handle invalid; skipped curve render.\n");
        }
    }
    else
    {
        printf("Curve JSON empty or invalid; skipped curve render: %s\n", curveJsonPath.c_str());
    }
    hostArena.Clear();
    device.deviceArena->Clear();

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
    return 0;
}
