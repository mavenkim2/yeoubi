#include "device/cuda_device.h"
#include "io/usd/load.h"
#include "scene/scene.h"
#include "tessellation/subdivision.h"
#include "texture/exr_io.h"
#include "texture/path_utils.h"
#include "texture/udim_utils.h"
#define STB_IMAGE_IMPLEMENTATION
#include "third_party/embree/tutorials/common/image/stb_image.h"
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
#include <cstring>
#include <filesystem>
#include <fstream>
#include <limits>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

#include <cuda_runtime_api.h>
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
    int materialTextureRefStride;
    int materialTextureRefSemanticCount;
    int textureViewSemantic;
    unsigned long long feedbackKeys;
    unsigned long long feedbackStats;
    int feedbackCapacity;
    int feedbackSamplePercent;
    int feedbackTileSize;
    int currentSpp;
};

enum class IntegratorType
{
    Primary,
    AO
};

enum class MaterialTextureSemantic : int
{
    Diffuse = 0,
    Roughness = 1,
    Metallic = 2,
    Occlusion = 3,
    Normal = 4,
    Ior = 5,
    Emissive = 6,
    Opacity = 7,
    Count = 8
};

struct CliOptions
{
    IntegratorType integrator = IntegratorType::Primary;
    MaterialTextureSemantic textureView = MaterialTextureSemantic::Diffuse;
    std::string inputPath;
    std::string outputPath = "optix_usd_scene.png";
    std::optional<ybi::float3> cameraPosition;
    std::optional<ybi::float3> lookAt;
    int spp = 1;
    bool useNtc = false;
    std::vector<std::string> purposes = {"default", "render"};
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

static bool TessellateRootSubdivisionMeshes(Scene *rootScene, const Camera &camera)
{
    YBI_ASSERT(rootScene);
    if (rootScene->subdivisionMeshes.empty())
    {
        return true;
    }
    if (!camera.hasValidCamera)
    {
        fprintf(stderr, "Subdivision requires a valid USD camera.\n");
        return false;
    }

    SubdivisionRunOptions options = {};
    options.level = 3;
    options.useCameraMatrices = true;
    options.cameraFromWorld = camera.cameraFromWorld;
    options.clipFromCamera = camera.clipFromCamera;
    options.viewportWidth = std::max(1, camera.viewportWidth);
    options.viewportHeight = std::max(1, camera.viewportHeight);
    options.verticalFovDegrees = camera.verticalFovDegrees;

    const size_t srcCount = rootScene->subdivisionMeshes.size();
    rootScene->meshes.reserve(rootScene->meshes.size() + srcCount);
    for (size_t i = 0; i < srcCount; ++i)
    {
        const SubdivisionMesh &subdivMesh = rootScene->subdivisionMeshes[i];
        SubdivisionRunResult result = {};
        if (!SubdivideAdaptive(subdivMesh, options, &result))
        {
            fprintf(stderr,
                    "Subdivision failed for rootScene->subdivisionMeshes[%zu] prim=%s\n",
                    i,
                    subdivMesh.primPath.c_str());
            return false;
        }

        result.mesh.materialIndex = subdivMesh.materialIndex;
        result.mesh.parentFromLocal = subdivMesh.parentFromLocal;
        rootScene->meshes.push_back(std::move(result.mesh));
    }
    rootScene->subdivisionMeshes.clear();
    return true;
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

static bool ParsePurposeToken(const std::string &input, std::string *outPurpose)
{
    YBI_ASSERT(outPurpose);
    *outPurpose = input;
    for (char &c : *outPurpose)
    {
        c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    }
    return *outPurpose == "default" || *outPurpose == "render" || *outPurpose == "proxy" ||
           *outPurpose == "guide";
}

static bool ParsePurposeList(const std::string &csv, std::vector<std::string> *outPurposes)
{
    YBI_ASSERT(outPurposes);
    outPurposes->clear();
    size_t begin = 0;
    while (begin <= csv.size())
    {
        const size_t comma = csv.find(',', begin);
        const size_t end = (comma == std::string::npos) ? csv.size() : comma;
        size_t b = begin;
        while (b < end && std::isspace(static_cast<unsigned char>(csv[b])))
        {
            ++b;
        }
        size_t e = end;
        while (e > b && std::isspace(static_cast<unsigned char>(csv[e - 1])))
        {
            --e;
        }
        std::string token = csv.substr(b, e - b);
        if (!token.empty())
        {
            std::string normalized;
            if (!ParsePurposeToken(token, &normalized))
            {
                return false;
            }
            if (std::find(outPurposes->begin(), outPurposes->end(), normalized) == outPurposes->end())
            {
                outPurposes->push_back(std::move(normalized));
            }
        }
        if (comma == std::string::npos)
        {
            break;
        }
        begin = comma + 1;
    }
    return !outPurposes->empty();
}

static std::string ToLowerCopy(const std::string &value)
{
    std::string out = value;
    for (char &c : out)
    {
        c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    }
    return out;
}

static bool ParseTextureViewSemantic(const std::string &input, MaterialTextureSemantic *outSemantic)
{
    YBI_ASSERT(outSemantic);
    const std::string token = ToLowerCopy(input);
    if (token == "diffuse")
    {
        *outSemantic = MaterialTextureSemantic::Diffuse;
        return true;
    }
    if (token == "roughness")
    {
        *outSemantic = MaterialTextureSemantic::Roughness;
        return true;
    }
    if (token == "metallic" || token == "metalness")
    {
        *outSemantic = MaterialTextureSemantic::Metallic;
        return true;
    }
    if (token == "occlusion" || token == "ao")
    {
        *outSemantic = MaterialTextureSemantic::Occlusion;
        return true;
    }
    if (token == "normal")
    {
        *outSemantic = MaterialTextureSemantic::Normal;
        return true;
    }
    if (token == "ior")
    {
        *outSemantic = MaterialTextureSemantic::Ior;
        return true;
    }
    if (token == "emissive" || token == "emissivecolor")
    {
        *outSemantic = MaterialTextureSemantic::Emissive;
        return true;
    }
    if (token == "opacity")
    {
        *outSemantic = MaterialTextureSemantic::Opacity;
        return true;
    }
    return false;
}

static bool TryMapInputNameToSemantic(const std::string &inputName, MaterialTextureSemantic *outSemantic)
{
    YBI_ASSERT(outSemantic);
    if (inputName == "diffuseColor")
    {
        *outSemantic = MaterialTextureSemantic::Diffuse;
        return true;
    }
    if (inputName == "roughness")
    {
        *outSemantic = MaterialTextureSemantic::Roughness;
        return true;
    }
    if (inputName == "metallic")
    {
        *outSemantic = MaterialTextureSemantic::Metallic;
        return true;
    }
    if (inputName == "occlusion")
    {
        *outSemantic = MaterialTextureSemantic::Occlusion;
        return true;
    }
    if (inputName == "normal")
    {
        *outSemantic = MaterialTextureSemantic::Normal;
        return true;
    }
    if (inputName == "ior")
    {
        *outSemantic = MaterialTextureSemantic::Ior;
        return true;
    }
    if (inputName == "emissiveColor")
    {
        *outSemantic = MaterialTextureSemantic::Emissive;
        return true;
    }
    if (inputName == "opacity")
    {
        *outSemantic = MaterialTextureSemantic::Opacity;
        return true;
    }
    return false;
}

static void PrintUsage(const char *exeName)
{
    printf("Usage: %s [--file path] [--out path] "
           "[--integrator primary|ao] [--spp N] "
           "[--cam-pos x y z] [--look-at x y z] [--ntc] [--view name] "
           "[--purposes csv] [--purpose name]\n",
           exeName);
    printf("  --file USDA/USD path\n");
    printf("  --out PNG output path\n");
    printf("  --integrator primary|ao\n");
    printf("  --spp spp passes; feedback dumped after each pass\n");
    printf("  --cam-pos optional camera position override\n");
    printf("  --look-at optional look-at target (default bounds center)\n");
    printf("  --ntc enable USD NTC decode path (falls back to image textures)\n");
    printf("  --view diffuse|roughness|metallic|occlusion|normal|ior|emissive|opacity\n");
    printf("  --purposes comma-separated default,render,proxy,guide\n");
    printf("  --purpose single purpose token, repeatable; overrides defaults\n");
}

static CliOptions ParseCli(int argc, char **argv)
{
    CliOptions options = {};
    bool purposeOverrideSet = false;
    for (int i = 1; i < argc; i++)
    {
        const std::string arg = argv[i];
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
        if (arg == "--ntc")
        {
            options.useNtc = true;
            continue;
        }
        if (arg == "--view")
        {
            if (i + 1 >= argc)
            {
                PrintUsage(argv[0]);
                std::abort();
            }
            MaterialTextureSemantic semantic = MaterialTextureSemantic::Diffuse;
            if (!ParseTextureViewSemantic(argv[++i], &semantic))
            {
                PrintUsage(argv[0]);
                std::abort();
            }
            options.textureView = semantic;
            continue;
        }
        if (arg == "--purposes")
        {
            if (i + 1 >= argc)
            {
                PrintUsage(argv[0]);
                std::abort();
            }
            std::vector<std::string> parsedPurposes;
            if (!ParsePurposeList(argv[++i], &parsedPurposes))
            {
                PrintUsage(argv[0]);
                std::abort();
            }
            options.purposes = std::move(parsedPurposes);
            purposeOverrideSet = true;
            continue;
        }
        if (arg == "--purpose")
        {
            if (i + 1 >= argc)
            {
                PrintUsage(argv[0]);
                std::abort();
            }
            std::string purpose;
            if (!ParsePurposeToken(argv[++i], &purpose))
            {
                PrintUsage(argv[0]);
                std::abort();
            }
            if (!purposeOverrideSet)
            {
                options.purposes.clear();
                purposeOverrideSet = true;
            }
            if (std::find(options.purposes.begin(), options.purposes.end(), purpose) ==
                options.purposes.end())
            {
                options.purposes.push_back(std::move(purpose));
            }
            continue;
        }

        PrintUsage(argv[0]);
        std::abort();
    }

    if (options.inputPath.empty())
    {
        PrintUsage(argv[0]);
        std::abort();
    }

    return options;
}

static bool SavePNG(const char *filePath, const std::vector<uint8_t> &rgba, int width, int height)
{
    const int strideInBytes = width * 4;
    return stbi_write_png(filePath, width, height, 4, rgba.data(), strideInBytes) != 0;
}

struct DecodedMaterialTexture
{
    bool valid = false;
    uint32_t udim = 1001u;
    int width = 0;
    int height = 0;
    std::vector<unsigned char> rgba8;
    std::string sourcePath;
    TextureWrapMode wrapS = TEXTURE_WRAP_MODE_REPEAT;
    TextureWrapMode wrapT = TEXTURE_WRAP_MODE_REPEAT;
};

struct UploadedMaterialTextures
{
    std::vector<LaunchParams::MaterialTextureRef> refs;
    std::vector<cudaArray_t> arrays;
    std::vector<cudaTextureObject_t> textureObjects;
};

static constexpr uint32_t kUdimMin = 1001u;
static constexpr uint32_t kUdimMax = 1100u;
static constexpr int kUdimSlotCount = 128;
static constexpr int kMaterialSemanticCount = static_cast<int>(MaterialTextureSemantic::Count);

static int MaterialTextureSlotIndex(
    size_t materialIndex, MaterialTextureSemantic semantic, uint32_t udim)
{
    int udimSlot = 0;
    if (udim >= kUdimMin && udim <= kUdimMax)
    {
        udimSlot = static_cast<int>(udim - kUdimMin);
    }
    udimSlot = std::max(0, std::min(udimSlot, kUdimSlotCount - 1));
    const int semanticIndex = static_cast<int>(semantic);
    return (static_cast<int>(materialIndex) * kMaterialSemanticCount + semanticIndex) *
               kUdimSlotCount +
           udimSlot;
}

static bool CheckCudaRuntime(cudaError_t result, std::string *outError, const char *callName)
{
    if (result == cudaSuccess)
    {
        return true;
    }
    if (outError)
    {
        *outError = std::string(callName) + " failed: " + cudaGetErrorString(result);
    }
    return false;
}

static cudaTextureAddressMode ToCudaAddressMode(TextureWrapMode mode)
{
    switch (mode)
    {
        case TEXTURE_WRAP_MODE_CLAMP:
            return cudaAddressModeClamp;
        case TEXTURE_WRAP_MODE_MIRROR:
            return cudaAddressModeMirror;
        case TEXTURE_WRAP_MODE_BLACK:
            return cudaAddressModeBorder;
        case TEXTURE_WRAP_MODE_REPEAT:
            return cudaAddressModeWrap;
        case TEXTURE_WRAP_MODE_USE_METADATA:
            return cudaAddressModeWrap;
        case TEXTURE_WRAP_MODE_UNKNOWN:
            return cudaAddressModeWrap;
    }
    return cudaAddressModeWrap;
}

static const MaterialTextureInput *
FindTextureInputBySemantic(const MaterialInfo &material, MaterialTextureSemantic semantic)
{
    for (const MaterialTextureInput &input : material.textureInputs)
    {
        if (input.texturePath.empty())
        {
            continue;
        }
        MaterialTextureSemantic inputSemantic = MaterialTextureSemantic::Diffuse;
        if (TryMapInputNameToSemantic(input.inputName, &inputSemantic) && inputSemantic == semantic)
        {
            return &input;
        }
    }
    return nullptr;
}

#if defined(YBI_OPTIX_HARNESS_WITH_NTC)
static bool LoadExrRgba8(const std::string &path,
                         int *outWidth,
                         int *outHeight,
                         std::vector<unsigned char> *outRgba8,
                         std::string *outReason)
{
    YBI_ASSERT(outWidth);
    YBI_ASSERT(outHeight);
    YBI_ASSERT(outRgba8);

    *outWidth = 0;
    *outHeight = 0;
    outRgba8->clear();

    std::vector<float> rgba;
    if (!ybi::texture::LoadExrRgba(path, outWidth, outHeight, &rgba, outReason, true))
    {
        return false;
    }

    const size_t pixelCount = static_cast<size_t>(*outWidth) * static_cast<size_t>(*outHeight);
    outRgba8->resize(pixelCount * 4u);
    for (size_t i = 0; i < pixelCount * 4u; ++i)
    {
        const float value = std::min(1.0f, std::max(0.0f, rgba[i]));
        (*outRgba8)[i] = static_cast<unsigned char>(value * 255.0f + 0.5f);
    }
    return true;
}
#endif

static bool LoadImageRgba8(const std::string &path,
                           int *outWidth,
                           int *outHeight,
                           std::vector<unsigned char> *outRgba8,
                           std::string *outReason)
{
    YBI_ASSERT(outWidth);
    YBI_ASSERT(outHeight);
    YBI_ASSERT(outRgba8);

    *outWidth = 0;
    *outHeight = 0;
    outRgba8->clear();

    if (ybi::texture::LowerExt(path) == ".exr")
    {
#if defined(YBI_OPTIX_HARNESS_WITH_NTC)
        return LoadExrRgba8(path, outWidth, outHeight, outRgba8, outReason);
#else
        if (outReason)
        {
            *outReason = "EXR loading unavailable (harness built without NTC/tinyexr): " + path;
        }
        return false;
#endif
    }

    int channels = 0;
    stbi_set_flip_vertically_on_load(1);
    stbi_uc *pixels = stbi_load(path.c_str(), outWidth, outHeight, &channels, 4);
    stbi_set_flip_vertically_on_load(0);
    if (!pixels || *outWidth <= 0 || *outHeight <= 0)
    {
        if (outReason)
        {
            *outReason = stbi_failure_reason() ? stbi_failure_reason() : "stbi_load failed";
        }
        if (pixels)
        {
            stbi_image_free(pixels);
        }
        return false;
    }

    const size_t bytes = static_cast<size_t>(*outWidth) * static_cast<size_t>(*outHeight) * 4u;
    outRgba8->resize(bytes);
    std::memcpy(outRgba8->data(), pixels, bytes);
    stbi_image_free(pixels);
    return true;
}

static bool DecodeImageTextures(const std::vector<MaterialInfo> &materials,
                                std::vector<DecodedMaterialTexture> *outTextures)
{
    YBI_ASSERT(outTextures);
    outTextures->clear();
    outTextures->resize(materials.size() * static_cast<size_t>(kMaterialSemanticCount) *
                        static_cast<size_t>(kUdimSlotCount));

    std::unordered_map<std::string, DecodedMaterialTexture> decodedByPath;
    int decodedTiles = 0;
    for (size_t materialIndex = 0; materialIndex < materials.size(); ++materialIndex)
    {
        const MaterialInfo &material = materials[materialIndex];
        for (const MaterialTextureInput &input : material.textureInputs)
        {
            if (input.texturePath.empty())
            {
                continue;
            }
            MaterialTextureSemantic semantic = MaterialTextureSemantic::Diffuse;
            if (!TryMapInputNameToSemantic(input.inputName, &semantic))
            {
                continue;
            }

            std::unordered_map<uint32_t, std::string> udimPaths;
            std::string udimReason;
            if (!ybi::usd_ntc::CollectUdimPaths(input.texturePath, udimPaths, udimReason))
            {
                std::printf("Image runtime: failed to resolve UDIMs for material %zu input %s (%s): %s\n",
                            materialIndex,
                            input.inputName.c_str(),
                            input.texturePath.c_str(),
                            udimReason.c_str());
                continue;
            }

            for (const auto &entry : udimPaths)
            {
                const uint32_t udim = entry.first;
                const std::string &tilePath = entry.second;
                const int dstIndex = MaterialTextureSlotIndex(materialIndex, semantic, udim);
                YBI_ASSERT(dstIndex >= 0 && static_cast<size_t>(dstIndex) < outTextures->size());

                auto cached = decodedByPath.find(tilePath);
                if (cached != decodedByPath.end())
                {
                    DecodedMaterialTexture tile = cached->second;
                    tile.udim = udim;
                    tile.wrapS = input.wrapS;
                    tile.wrapT = input.wrapT;
                    (*outTextures)[static_cast<size_t>(dstIndex)] = std::move(tile);
                    decodedTiles++;
                    continue;
                }

                DecodedMaterialTexture texture = {};
                texture.valid = true;
                texture.udim = udim;
                texture.sourcePath = tilePath;
                texture.wrapS = input.wrapS;
                texture.wrapT = input.wrapT;

                std::string reason;
                if (!LoadImageRgba8(tilePath, &texture.width, &texture.height, &texture.rgba8, &reason))
                {
                    std::printf("Image runtime: failed to load material %zu input %s tile %u (%s): %s\n",
                                materialIndex,
                                input.inputName.c_str(),
                                udim,
                                tilePath.c_str(),
                                reason.c_str());
                    continue;
                }

                (*outTextures)[static_cast<size_t>(dstIndex)] = texture;
                decodedByPath.emplace(tilePath, std::move(texture));
                decodedTiles++;
            }
        }
    }

    std::printf("Image runtime: decoded %d/%zu material UDIM tiles\n",
                decodedTiles,
                materials.size() * static_cast<size_t>(kMaterialSemanticCount) *
                    static_cast<size_t>(kUdimSlotCount));
    return true;
}

static bool UploadDecodedTexturesToCuda(const std::vector<DecodedMaterialTexture> &decodedTextures,
                                        UploadedMaterialTextures *outTextures,
                                        std::string *outError)
{
    YBI_ASSERT(outTextures);
    for (cudaTextureObject_t textureObject : outTextures->textureObjects)
    {
        if (textureObject)
        {
            cudaDestroyTextureObject(textureObject);
        }
    }
    for (cudaArray_t array : outTextures->arrays)
    {
        if (array)
        {
            cudaFreeArray(array);
        }
    }
    outTextures->arrays.clear();
    outTextures->textureObjects.clear();
    outTextures->refs.clear();
    outTextures->refs.resize(decodedTextures.size());

    for (size_t i = 0; i < decodedTextures.size(); ++i)
    {
        const DecodedMaterialTexture &src = decodedTextures[i];
        if (!src.valid || src.width <= 0 || src.height <= 0 || src.rgba8.empty())
        {
            continue;
        }

        cudaChannelFormatDesc channelDesc =
            cudaCreateChannelDesc(8, 8, 8, 8, cudaChannelFormatKindUnsigned);
        cudaArray_t array = nullptr;
        if (!CheckCudaRuntime(
                cudaMallocArray(&array, &channelDesc, src.width, src.height), outError, "cudaMallocArray"))
        {
            return false;
        }

        if (!CheckCudaRuntime(cudaMemcpy2DToArray(array,
                                                  0,
                                                  0,
                                                  src.rgba8.data(),
                                                  static_cast<size_t>(src.width) * 4u,
                                                  static_cast<size_t>(src.width) * 4u,
                                                  src.height,
                                                  cudaMemcpyHostToDevice),
                              outError,
                              "cudaMemcpy2DToArray"))
        {
            cudaFreeArray(array);
            return false;
        }

        cudaResourceDesc resourceDesc = {};
        resourceDesc.resType = cudaResourceTypeArray;
        resourceDesc.res.array.array = array;

        cudaTextureDesc textureDesc = {};
        textureDesc.addressMode[0] = ToCudaAddressMode(src.wrapS);
        textureDesc.addressMode[1] = ToCudaAddressMode(src.wrapT);
        textureDesc.filterMode = cudaFilterModePoint;
        textureDesc.readMode = cudaReadModeNormalizedFloat;
        textureDesc.normalizedCoords = 1;

        cudaTextureObject_t textureObject = 0;
        if (!CheckCudaRuntime(cudaCreateTextureObject(&textureObject, &resourceDesc, &textureDesc, nullptr),
                              outError,
                              "cudaCreateTextureObject"))
        {
            cudaFreeArray(array);
            return false;
        }

        outTextures->arrays.push_back(array);
        outTextures->textureObjects.push_back(textureObject);
        outTextures->refs[i].textureObject = static_cast<unsigned long long>(textureObject);
        outTextures->refs[i].width = src.width;
        outTextures->refs[i].height = src.height;
        outTextures->refs[i].valid = 1;
    }

    return true;
}

static void DestroyUploadedTextures(UploadedMaterialTextures *textures)
{
    if (!textures)
    {
        return;
    }
    for (cudaTextureObject_t textureObject : textures->textureObjects)
    {
        if (textureObject)
        {
            cudaDestroyTextureObject(textureObject);
        }
    }
    for (cudaArray_t array : textures->arrays)
    {
        if (array)
        {
            cudaFreeArray(array);
        }
    }
    textures->textureObjects.clear();
    textures->arrays.clear();
    textures->refs.clear();
}

static unsigned int FeedbackTileX(unsigned long long key)
{
    return static_cast<unsigned int>((key >> 0u) & 0x1ffull);
}

static unsigned int FeedbackTileY(unsigned long long key)
{
    return static_cast<unsigned int>((key >> 9u) & 0x1ffull);
}

static unsigned int FeedbackUdim(unsigned long long key)
{
    return 1001u + static_cast<unsigned int>((key >> 18u) & 0x7full);
}

static unsigned int FeedbackTextureId(unsigned long long key)
{
    return static_cast<unsigned int>((key >> 25u) & 0x7fffffull);
}

static unsigned int FeedbackMip(unsigned long long key)
{
    return static_cast<unsigned int>((key >> 48u) & 0xfull);
}

struct UploadedMeshRefs
{
    std::vector<LaunchParams::InstanceGeomRef> refs;
    std::vector<DeviceMemoryView<uint8_t>> ownedBuffers;
};

static const Attribute *FindMeshSTAttribute(const Mesh &mesh)
{
    for (const Attribute &attr : mesh.attributes)
    {
        if (attr.name == "st" && attr.type == AttributeType::Float2 &&
            attr.interpolation == PrimvarInterpolation::FaceVarying && attr.indices.size() > 0)
        {
            if ((attr.data.size() % sizeof(ybi::float2)) == 0)
            {
                return &attr;
            }
        }
    }
    return nullptr;
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
        int texcoordCount = 0;
        int texcoordIndexCount = 0;
        const Attribute *stAttr = FindMeshSTAttribute(mesh);
        if (stAttr)
        {
            const size_t texcoordsBytes = stAttr->data.size();
            const size_t texcoordIndicesBytes = sizeof(int) * stAttr->indices.size();
            texcoordsBuffer = device->AllocBytes(texcoordsBytes);
            texcoordIndicesBuffer = device->AllocBytes(texcoordIndicesBytes);
            device->CopyBytesToDevice(texcoordsBuffer, stAttr->data.data(), texcoordsBytes);
            device->CopyBytesToDevice(texcoordIndicesBuffer, stAttr->indices.data(), texcoordIndicesBytes);
            out.ownedBuffers.push_back(texcoordsBuffer);
            out.ownedBuffers.push_back(texcoordIndicesBuffer);
            texcoordCount = int(stAttr->data.size() / sizeof(ybi::float2));
            texcoordIndexCount = int(stAttr->indices.size());
        }

        out.refs[job.refIndex] = {(unsigned long long)positionsBuffer.data(),
                                  (unsigned long long)indicesBuffer.data(),
                                  (unsigned long long)texcoordsBuffer.data(),
                                  (unsigned long long)texcoordIndicesBuffer.data(),
                                  (int)mesh.positions.size(),
                                  (int)mesh.indices.size(),
                                  texcoordCount,
                                  texcoordIndexCount,
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
                              int materialTextureRefStride,
                              int materialTextureRefSemanticCount,
                              MaterialTextureSemantic textureViewSemantic,
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

    const int sppPassCount = std::max(1, spp);
    const int feedbackCapacity = std::max(1, width * height);
    const size_t feedbackKeysBytes =
        static_cast<size_t>(feedbackCapacity) * sizeof(unsigned long long);
    const size_t feedbackStatsBytes = 2u * sizeof(unsigned int);
    CUdeviceptr feedbackKeysBuffer = 0;
    CUdeviceptr feedbackStatsBuffer = 0;
    CUDA_ASSERT(cuMemAlloc(&feedbackKeysBuffer, feedbackKeysBytes));
    CUDA_ASSERT(cuMemAlloc(&feedbackStatsBuffer, feedbackStatsBytes));

    std::filesystem::path feedbackDir = std::filesystem::path(outputFile);
    feedbackDir += ".feedback";
    std::error_code feedbackEc;
    std::filesystem::create_directories(feedbackDir, feedbackEc);
    if (feedbackEc)
    {
        fprintf(stderr, "Failed to create feedback dir: %s\n", feedbackDir.string().c_str());
        std::abort();
    }

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
    params.materialTextureRefStride = materialTextureRefStride;
    params.materialTextureRefSemanticCount = materialTextureRefSemanticCount;
    params.textureViewSemantic = static_cast<int>(textureViewSemantic);
    params.feedbackKeys = (unsigned long long)feedbackKeysBuffer;
    params.feedbackStats = (unsigned long long)feedbackStatsBuffer;
    params.feedbackCapacity = feedbackCapacity;
    params.feedbackSamplePercent = 10;
    params.feedbackTileSize = 128;
    params.currentSpp = 0;

    CUdeviceptr paramsBuffer = 0;
    CUDA_ASSERT(cuMemAlloc(&paramsBuffer, sizeof(LaunchParams)));
    printf("render: params buffer allocated\n");
    fflush(stdout);

    std::vector<uint8_t> hostPassImage(imageSize, 0);
    std::vector<float> accumRgb(static_cast<size_t>(width) * static_cast<size_t>(height) * 3u, 0.0f);
    std::vector<unsigned long long> feedbackKeysHost(feedbackCapacity, 0ull);
    unsigned int feedbackStatsHost[2] = {0u, 0u};
    unsigned int feedbackStatsZero[2] = {0u, 0u};

    for (int sppIndex = 0; sppIndex < sppPassCount; ++sppIndex)
    {
        params.spp = integrator == IntegratorType::AO ? 1 : 1;
        params.currentSpp = sppIndex;
        CUDA_ASSERT(cuMemcpyHtoD(paramsBuffer, &params, sizeof(LaunchParams)));
        CUDA_ASSERT(cuMemcpyHtoD(feedbackStatsBuffer, feedbackStatsZero, feedbackStatsBytes));

        OPTIX_CHECK(
            optixLaunch(pipeline, 0, paramsBuffer, sizeof(LaunchParams), &sbt, width, height, 1));
        CUDA_ASSERT(cuStreamSynchronize(0));

        CUDA_ASSERT(cuMemcpyDtoH(hostPassImage.data(), imageBuffer, imageSize));
        const size_t pixelCount = static_cast<size_t>(width) * static_cast<size_t>(height);
        for (size_t i = 0; i < pixelCount; ++i)
        {
            accumRgb[i * 3 + 0] += hostPassImage[i * 4 + 0] / 255.0f;
            accumRgb[i * 3 + 1] += hostPassImage[i * 4 + 1] / 255.0f;
            accumRgb[i * 3 + 2] += hostPassImage[i * 4 + 2] / 255.0f;
        }

        CUDA_ASSERT(cuMemcpyDtoH(feedbackStatsHost, feedbackStatsBuffer, feedbackStatsBytes));
        const unsigned int sampledCount = feedbackStatsHost[0];
        const unsigned int overflowCount = feedbackStatsHost[1];
        const unsigned int copyCount = std::min(sampledCount, (unsigned int)feedbackCapacity);
        if (copyCount > 0)
        {
            CUDA_ASSERT(cuMemcpyDtoH(
                feedbackKeysHost.data(), feedbackKeysBuffer, size_t(copyCount) * sizeof(unsigned long long)));
        }

        std::unordered_map<unsigned long long, unsigned int> histogram;
        histogram.reserve(copyCount);
        for (unsigned int i = 0; i < copyCount; ++i)
        {
            histogram[feedbackKeysHost[i]] += 1u;
        }

        char feedbackFileName[256];
        std::snprintf(feedbackFileName, sizeof(feedbackFileName), "spp_%04d.txt", sppIndex);
        const std::filesystem::path feedbackPath = feedbackDir / feedbackFileName;
        std::FILE *feedbackFile = std::fopen(feedbackPath.string().c_str(), "w");
        if (!feedbackFile)
        {
            fprintf(stderr, "Failed to write feedback file: %s\n", feedbackPath.string().c_str());
            std::abort();
        }
        std::fprintf(feedbackFile,
                     "spp=%d sampled=%u stored=%u overflow=%u unique=%zu\n",
                     sppIndex,
                     sampledCount,
                     copyCount,
                     overflowCount,
                     histogram.size());
        std::fprintf(feedbackFile, "textureId udim tileX tileY mip count\n");
        for (const auto &it : histogram)
        {
            const unsigned long long key = it.first;
            std::fprintf(feedbackFile,
                         "%u %u %u %u %u %u\n",
                         FeedbackTextureId(key),
                         FeedbackUdim(key),
                         FeedbackTileX(key),
                         FeedbackTileY(key),
                         FeedbackMip(key),
                         it.second);
        }
        std::fclose(feedbackFile);
    }

    std::vector<uint8_t> hostImage(imageSize, 0);
    const float invSpp = 1.0f / float(sppPassCount);
    for (size_t i = 0, pixelCount = static_cast<size_t>(width) * static_cast<size_t>(height); i < pixelCount;
         ++i)
    {
        const float r = std::min(1.0f, std::max(0.0f, accumRgb[i * 3 + 0] * invSpp));
        const float g = std::min(1.0f, std::max(0.0f, accumRgb[i * 3 + 1] * invSpp));
        const float b = std::min(1.0f, std::max(0.0f, accumRgb[i * 3 + 2] * invSpp));
        hostImage[i * 4 + 0] = static_cast<uint8_t>(r * 255.0f + 0.5f);
        hostImage[i * 4 + 1] = static_cast<uint8_t>(g * 255.0f + 0.5f);
        hostImage[i * 4 + 2] = static_cast<uint8_t>(b * 255.0f + 0.5f);
        hostImage[i * 4 + 3] = 255u;
    }

    const bool writeOk = SavePNG(outputFile, hostImage, width, height);
    printf("render: save returned=%d file=%s\n", writeOk ? 1 : 0, outputFile);
    fflush(stdout);
    if (!writeOk)
    {
        fprintf(stderr, "Failed to write PNG: %s\n", outputFile);
    }

    CUDA_ASSERT(cuMemFree(paramsBuffer));
    CUDA_ASSERT(cuMemFree(feedbackStatsBuffer));
    CUDA_ASSERT(cuMemFree(feedbackKeysBuffer));
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

    printf("optix_harness: loading usd scene %s\n", options.inputPath.c_str());
    fflush(stdout);

    ScenePool scenePool = {};
    USDLoadOptions loadOptions = {};
    loadOptions.purposes = options.purposes;
    LoadUSDScene(&scenePool, options.inputPath, loadOptions);
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
    Scene *rootScene = flattenedScenePool.scenes[flattenedScenePool.rootSceneIndex].get();
    YBI_ASSERT(rootScene);
    if (!TessellateRootSubdivisionMeshes(rootScene, flattenedScenePool.camera))
    {
        return 1;
    }

    DeviceMemoryView<uint8_t> materialTextureRefsBuffer = {};
    int materialTextureRefCount = 0;
    int materialTextureRefStride = kUdimSlotCount;
    int materialTextureRefSemanticCount = kMaterialSemanticCount;
    UploadedMaterialTextures uploadedMaterialTextures = {};
    std::vector<DecodedMaterialTexture> decodedTextures;
    if (!DecodeImageTextures(scenePool.materials, &decodedTextures))
    {
        fprintf(stderr, "Image runtime decode failed.\n");
        return 1;
    }

    if (options.useNtc)
    {
#if defined(YBI_OPTIX_HARNESS_WITH_NTC)
        std::vector<testbvh::DecodedMaterialTexture> ntcTextures;
        std::string ntcError;
        if (testbvh::DecodeNtcDiffuseTextures(scenePool.materials, &ntcTextures, &ntcError))
        {
            int overrideCount = 0;
            const size_t numMaterials = std::min(scenePool.materials.size(), ntcTextures.size());
            for (size_t i = 0; i < numMaterials; ++i)
            {
                if (!ntcTextures[i].valid)
                {
                    continue;
                }
                const int slot =
                    MaterialTextureSlotIndex(i, MaterialTextureSemantic::Diffuse, kUdimMin);
                YBI_ASSERT(slot >= 0 && static_cast<size_t>(slot) < decodedTextures.size());
                DecodedMaterialTexture &dst = decodedTextures[static_cast<size_t>(slot)];
                dst.valid = true;
                dst.udim = kUdimMin;
                dst.width = ntcTextures[i].width;
                dst.height = ntcTextures[i].height;
                dst.rgba8 = std::move(ntcTextures[i].rgba8);
                dst.sourcePath = ntcTextures[i].ntcPath;
                const MaterialTextureInput *diffuse =
                    FindTextureInputBySemantic(scenePool.materials[i], MaterialTextureSemantic::Diffuse);
                dst.wrapS = diffuse ? diffuse->wrapS : TEXTURE_WRAP_MODE_REPEAT;
                dst.wrapT = diffuse ? diffuse->wrapT : TEXTURE_WRAP_MODE_REPEAT;
                overrideCount++;
            }
            printf("NTC runtime: applied %d material overrides\n", overrideCount);
        }
        else
        {
            fprintf(stderr, "NTC runtime decode failed (continuing with image textures): %s\n", ntcError.c_str());
        }
#else
        fprintf(stderr,
                "NTC runtime requested via --ntc, but harness built without WITH_NTC. "
                "Using image textures only.\n");
#endif
    }

    std::string textureUploadError;
    if (!UploadDecodedTexturesToCuda(decodedTextures, &uploadedMaterialTextures, &textureUploadError))
    {
        fprintf(stderr, "Texture upload failed: %s\n", textureUploadError.c_str());
        return 1;
    }

    if (!uploadedMaterialTextures.refs.empty())
    {
        std::vector<LaunchParams::MaterialTextureRef> launchRefs(uploadedMaterialTextures.refs.size());
        for (size_t i = 0; i < uploadedMaterialTextures.refs.size(); ++i)
        {
            const LaunchParams::MaterialTextureRef &src = uploadedMaterialTextures.refs[i];
            launchRefs[i] = {src.textureObject, src.width, src.height, src.valid};
        }
        const size_t refsBytes = launchRefs.size() * sizeof(LaunchParams::MaterialTextureRef);
        materialTextureRefsBuffer = device.AllocBytes(refsBytes);
        device.CopyBytesToDevice(materialTextureRefsBuffer, launchRefs.data(), refsBytes);
        materialTextureRefCount = static_cast<int>(scenePool.materials.size());
    }

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
                      materialTextureRefStride,
                      materialTextureRefSemanticCount,
                      options.textureView,
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
    DestroyUploadedTextures(&uploadedMaterialTextures);
    for (DeviceMemoryView<uint8_t> &buffer : uploadedRefs.ownedBuffers)
    {
        device.FreeBytes(buffer);
    }
    hostArena.Clear();
    device.deviceArena->Clear();

    device.DestroyOptixPrimaryPipeline();
    std::fflush(stdout);
    std::fflush(stderr);
    std::_Exit(0);
}
