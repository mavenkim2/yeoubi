#include "device/cuda_device.h"
#include "io/usd/load.h"
#include "../io/curve_json_io.h"
#include "../io/obj_mesh_io.h"
#include "scene/scene.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "third_party/stb_image_write.h"
#include "util/array.h"
#include "util/float3.h"
#include "util/float3x4.h"
#include "util/float4.h"
#include "util/float4x4.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <optional>
#include <string>
#include <vector>

#include <optix_stubs.h>

using namespace ybi;

#if defined(YBI_OPTIX_HARNESS_WITH_NTC)
#include "optix_ntc_runtime.h"
#endif

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
        unsigned long long texcoords;
        unsigned long long texcoordIndices;
        int numPositions;
        int numIndices;
        int numTexcoords;
        int numTexcoordIndices;
        int materialIndex;
    };

    struct WireframeConfig
    {
        float lineWidth;
        float lineFeather;
        float edgeDarkness;
        float padding;
    };

    struct MaterialTextureRef
    {
        unsigned long long textureObject;
        int width;
        int height;
        int valid;
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
    unsigned long long materialTextureRefs;
    int materialTextureRefCount;
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
    std::string ntcDir;
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
           "[--integrator primary|ao] [--spp N] [--cam-pos x y z] [--look-at x y z] [--ntc-dir path]\n",
           exeName);
    printf("  --type triangle|cluster|curve|usd\n");
    printf("  --file OBJ path for triangle/cluster; JSON path for curve; USDA/USD path for usd\n");
    printf("  --out PNG output path\n");
    printf("  --integrator primary|ao\n");
    printf("  --spp samples per pixel for ao integrator\n");
    printf("  --cam-pos optional camera position override\n");
    printf("  --look-at optional look-at target (default bounds center)\n");
    printf("  --ntc-dir optional directory with <sanitizedMaterialPath>.ntc files\n");
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
        if (arg == "--ntc-dir")
        {
            if (i + 1 >= argc)
            {
                PrintUsage(argv[0]);
                std::abort();
            }
            options.ntcDir = argv[++i];
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

static constexpr unsigned int kHitgroupSbtOffset = 0u;

static ybi::float3x4 IdentityTransform3x4()
{
    return ybi::float3x4(1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f);
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

struct UploadedMeshRefs
{
    std::vector<LaunchParams::InstanceGeomRef> refs;
    std::vector<DeviceMemoryView<uint8_t>> ownedBuffers;
};

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

static UploadedMeshRefs UploadScenePoolMeshRefs(
    CUDADevice *device, const std::vector<SceneMeshUploadRef> &meshUploadRefs)
{
    YBI_ASSERT(device);
    UploadedMeshRefs out;
    out.refs.resize(meshUploadRefs.size());

    for (const SceneMeshUploadRef &job : meshUploadRefs)
    {
        YBI_ASSERT(job.mesh);
        YBI_ASSERT(job.refIndex < out.refs.size());
        const Mesh &mesh = *job.mesh;
        if (mesh.positions.size() == 0 || mesh.indices.size() == 0)
        {
            continue;
        }

        const size_t positionsBytes = sizeof(ybi::float3) * mesh.positions.size();
        const size_t indicesBytes = sizeof(int) * mesh.indices.size();
        DeviceMemoryView<uint8_t> positionsBuffer = device->AllocBytes(positionsBytes);
        DeviceMemoryView<uint8_t> indicesBuffer = device->AllocBytes(indicesBytes);
        device->CopyBytesToDevice(positionsBuffer, mesh.positions.data(), positionsBytes);
        device->CopyBytesToDevice(indicesBuffer, mesh.indices.data(), indicesBytes);

        DeviceMemoryView<uint8_t> texcoordsBuffer = {};
        DeviceMemoryView<uint8_t> texcoordIndicesBuffer = {};
        if (mesh.texcoords.size() > 0 && mesh.texcoordIndices.size() > 0)
        {
            const size_t texcoordsBytes = sizeof(ybi::float2) * mesh.texcoords.size();
            const size_t texcoordIndicesBytes = sizeof(int) * mesh.texcoordIndices.size();
            texcoordsBuffer = device->AllocBytes(texcoordsBytes);
            texcoordIndicesBuffer = device->AllocBytes(texcoordIndicesBytes);
            device->CopyBytesToDevice(texcoordsBuffer, mesh.texcoords.data(), texcoordsBytes);
            device->CopyBytesToDevice(
                texcoordIndicesBuffer, mesh.texcoordIndices.data(), texcoordIndicesBytes);
            out.ownedBuffers.push_back(texcoordsBuffer);
            out.ownedBuffers.push_back(texcoordIndicesBuffer);
        }

        out.refs[job.refIndex] = {(unsigned long long)positionsBuffer.data(),
                                  (unsigned long long)indicesBuffer.data(),
                                  (unsigned long long)texcoordsBuffer.data(),
                                  (unsigned long long)texcoordIndicesBuffer.data(),
                                  (int)mesh.positions.size(),
                                  (int)mesh.indices.size(),
                                  (int)mesh.texcoords.size(),
                                  (int)mesh.texcoordIndices.size(),
                                  mesh.materialIndex};
        out.ownedBuffers.push_back(positionsBuffer);
        out.ownedBuffers.push_back(indicesBuffer);
    }
    return out;
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
                              CUdeviceptr materialTextureRefsBuffer,
                              int materialTextureRefCount,
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
    params.materialTextureRefs = (unsigned long long)materialTextureRefsBuffer;
    params.materialTextureRefCount = materialTextureRefCount;

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
        Mesh mesh = testio::LoadObjMesh(options.inputPath);
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
            0ull,
            0ull,
            (int)mesh.positions.size(),
            (int)mesh.indices.size(),
            0,
            0,
            -1,
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
                                  0,
                                  0,
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
            if (!device.SupportsClusterAccel())
            {
                printf("Cluster GAS not supported on this OptiX device/context.\n");
            }
            else
            {
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
                                      0,
                                      0,
                                      std::nullopt,
                                      options.cameraPosition,
                                      options.lookAt);
                    printf("Wrote %s\n", options.outputPath.c_str());
                }
                else
                {
                    printf("Cluster handle invalid; skipped cluster render.\n");
                }
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
        Curves curves = testio::LoadCurveJson(options.inputPath);
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

        ScenePool flattenedScenePool = {};
        std::string flattenError;
        if (!FlattenScenePoolToRootChildren(&scenePool, &flattenedScenePool, &flattenError))
        {
            fprintf(stderr, "Failed to flatten USD ScenePool: %s\n", flattenError.c_str());
            return 1;
        }
        if (flattenedScenePool.scenes.empty() ||
            flattenedScenePool.rootSceneIndex >= flattenedScenePool.scenes.size())
        {
            fprintf(stderr,
                    "Flattened ScenePool invalid for USD scene: %s\n",
                    options.inputPath.c_str());
            return 1;
        }

        DeviceMemoryView<uint8_t> materialTextureRefsBuffer = {};
        int materialTextureRefCount = 0;
#if defined(YBI_OPTIX_HARNESS_WITH_NTC)
        testbvh::UploadedMaterialTextures uploadedMaterialTextures = {};
        if (!options.ntcDir.empty())
        {
            std::vector<testbvh::DecodedMaterialTexture> decodedTextures;
            std::string ntcError;
            const bool ok = testbvh::DecodeNtcDiffuseTextures(
                scenePool.materials, options.ntcDir, &decodedTextures, &ntcError);
            if (!ok)
            {
                fprintf(stderr, "NTC runtime decode failed: %s\n", ntcError.c_str());
                return 1;
            }

            if (!testbvh::UploadDecodedTexturesToCuda(
                    decodedTextures, &uploadedMaterialTextures, &ntcError))
            {
                fprintf(stderr, "NTC runtime upload failed: %s\n", ntcError.c_str());
                return 1;
            }
            if (!uploadedMaterialTextures.refs.empty())
            {
                std::vector<LaunchParams::MaterialTextureRef> launchRefs(
                    uploadedMaterialTextures.refs.size());
                for (size_t i = 0; i < uploadedMaterialTextures.refs.size(); ++i)
                {
                    const testbvh::MaterialTextureRef &src = uploadedMaterialTextures.refs[i];
                    launchRefs[i] = {src.textureObject, src.width, src.height, src.valid};
                }
                const size_t refsBytes =
                    launchRefs.size() * sizeof(LaunchParams::MaterialTextureRef);
                materialTextureRefsBuffer = device.AllocBytes(refsBytes);
                device.CopyBytesToDevice(materialTextureRefsBuffer, launchRefs.data(), refsBytes);
                materialTextureRefCount = static_cast<int>(launchRefs.size());
            }
        }
#else
        if (!options.ntcDir.empty())
        {
            fprintf(stderr,
                    "NTC runtime decode requested but harness was built without WITH_NTC support.\n");
            return 1;
        }
#endif

        Scene *rootScene = flattenedScenePool.scenes[flattenedScenePool.rootSceneIndex].get();
        YBI_ASSERT(rootScene);

        std::vector<SceneMeshUploadRef> meshUploadRefs;
        CollectScenePoolMeshUploadRefs(&flattenedScenePool, &meshUploadRefs);

        for (const std::unique_ptr<Scene> &scenePtr : flattenedScenePool.scenes)
        {
            Scene *scene = scenePtr.get();
            YBI_ASSERT(scene);
            device.BuildBVH(scene);
        }
        if (rootScene->bvhHandle == 0)
        {
            fprintf(stderr, "USD root scene BVH is invalid: %s\n", options.inputPath.c_str());
            return 1;
        }

        UploadedMeshRefs uploadedRefs = UploadScenePoolMeshRefs(&device, meshUploadRefs);

        DeviceMemoryView<uint8_t> instanceGeomRefsBuffer = {};
        if (!uploadedRefs.refs.empty())
        {
            const size_t refsBytes =
                uploadedRefs.refs.size() * sizeof(LaunchParams::InstanceGeomRef);
            instanceGeomRefsBuffer = device.AllocBytes(refsBytes);
            device.CopyBytesToDevice(instanceGeomRefsBuffer, uploadedRefs.refs.data(), refsBytes);
        }

        std::optional<RenderCameraOverride> usdCamera = std::nullopt;
        if (!options.cameraPosition.has_value() && !options.lookAt.has_value())
        {
            usdCamera = BuildUsdRenderCamera(flattenedScenePool.camera);
            if (usdCamera.has_value())
            {
                printf("optix_harness: using usd camera viewport=%dx%d\n",
                       usdCamera->width,
                       usdCamera->height);
            }
            else
            {
                fprintf(stderr, "USD camera missing/invalid: %s\n", options.inputPath.c_str());
                return 1;
            }
            fflush(stdout);
        }

        const ybi::float3 dummyBoundsMin = ybi::make_float3(-1.0f, -1.0f, -1.0f);
        const ybi::float3 dummyBoundsMax = ybi::make_float3(1.0f, 1.0f, 1.0f);

        RenderTraversable(pipeline,
                          sbt,
                          (OptixTraversableHandle)rootScene->bvhHandle,
                          dummyBoundsMin,
                          dummyBoundsMax,
                          options.outputPath.c_str(),
                          options.integrator,
                          options.spp,
                          (CUdeviceptr)instanceGeomRefsBuffer.data(),
                          (int)uploadedRefs.refs.size(),
                          (CUdeviceptr)materialTextureRefsBuffer.data(),
                          materialTextureRefCount,
                          usdCamera,
                          options.cameraPosition,
                          options.lookAt);
        printf("Wrote %s\n", options.outputPath.c_str());
        if (instanceGeomRefsBuffer.data())
        {
            device.FreeBytes(instanceGeomRefsBuffer);
        }
        if (materialTextureRefsBuffer.data())
        {
            device.FreeBytes(materialTextureRefsBuffer);
        }
#if defined(YBI_OPTIX_HARNESS_WITH_NTC)
        testbvh::DestroyUploadedTextures(&uploadedMaterialTextures);
#endif
        for (DeviceMemoryView<uint8_t> &buffer : uploadedRefs.ownedBuffers)
        {
            device.FreeBytes(buffer);
        }
    }
    hostArena.Clear();
    device.deviceArena->Clear();

    device.DestroyOptixPrimaryPipeline();
    std::fflush(stdout);
    std::fflush(stderr);
    std::_Exit(0);
}
