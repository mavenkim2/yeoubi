#include "device/device.h"
#include "io/usd/load.h"
#include "render/dispatch_types.h"
#include "render/integrator_ray_differential.h"
#include "render/launch_params.h"
#include "render/scene_memory.h"
#include "scene/scene.h"
#include "tessellation/subdivision.h"
#include "texture/exr_io.h"
#include "texture/path_utils.h"
#include "texture/udim_utils.h"
#include "texture/virtual_texture/feedback.h"
#include "texture/virtual_texture/helpers.h"
#include "texture/virtual_texture/manager.h"
#include "texture/virtual_texture/material_metadata.h"
#include "texture/virtual_texture/tile_file.h"
#define STB_IMAGE_IMPLEMENTATION
#include "third_party/embree/tutorials/common/image/stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "third_party/stb_image_write.h"
#include "util/array.h"
#include "util/half_float.h"
#include "util/vec3.h"
#include "util/float3x4.h"
#include "util/vec4.h"
#include "util/float4x4.h"
#include <algorithm>
#include <atomic>
#include <cctype>
#include <chrono>
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
#include <unordered_set>
#include <vector>
#include <tbb/parallel_for.h>

#if defined(YBI_OPTIX_HARNESS_WITH_NTC)
#include "../../tests/bvh/optix/optix_ntc_runtime.h"
#endif

using namespace ybi;

#ifndef YBI_OPTIX_PRIMARY_PTX_PATH
#define YBI_OPTIX_PRIMARY_PTX_PATH ""
#endif

namespace
{
enum class IntegratorType
{
    Primary,
    AO,
    Path
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
    SpecularColor = 8,
    Clearcoat = 9,
    ClearcoatRoughness = 10,
    Count = 11
};

struct CliOptions
{
    IntegratorType integrator = IntegratorType::Primary;
    MaterialTextureSemantic textureView = MaterialTextureSemantic::Diffuse;
    DeviceKind deviceKind = DeviceKind::GPU;
    std::string inputPath;
    std::string outputPath = "optix_usd_scene.png";
    std::optional<ybi::Vec3> cameraPosition;
    std::optional<ybi::Vec3> lookAt;
    std::string usdCamera;
    int spp = 1;
    int maxDepth = 4;
    bool useNtc = false;
    bool virtualTexture = false;
    std::string virtualTextureTilesDir;
    uint64_t virtualTextureCacheBytes = 0u;
    int virtualTextureTailMaxDim = 32;
    int virtualTextureMaxPageUploads = 0;
    bool singlePixelMode = false;
    int singlePixelX = 0;
    int singlePixelY = 0;
    bool writeFeedbackFiles = false;
    std::vector<std::string> purposes = {"default", "render"};
};

struct RenderCameraOverride
{
    int width = 1280;
    int height = 720;
    ybi::Vec3 origin = ybi::Vec3(0.0f, 0.0f, 0.0f);
    ybi::Float4x4 cameraFromWorld = ybi::Float4x4::Identity();
    ybi::Float4x4 worldFromCamera = ybi::Float4x4::Identity();
    ybi::Float4x4 clipFromCamera = ybi::Float4x4::Identity();
    ybi::Float4x4 cameraFromRaster = ybi::Float4x4::Identity();
    ybi::Float4x4 rasterFromCamera = ybi::Float4x4::Identity();
    ybi::Vec3 minPosDifferentialX = ybi::Vec3(0.0f);
    ybi::Vec3 minPosDifferentialY = ybi::Vec3(0.0f);
    ybi::Vec3 minDirDifferentialX = ybi::Vec3(0.0f);
    ybi::Vec3 minDirDifferentialY = ybi::Vec3(0.0f);
    ybi::Vec3 U = ybi::Vec3(1.0f, 0.0f, 0.0f);
    ybi::Vec3 V = ybi::Vec3(0.0f, 1.0f, 0.0f);
    ybi::Vec3 W = ybi::Vec3(0.0f, 0.0f, 1.0f);
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

static bool InvertAffine(const ybi::Float4x4 &m, ybi::Float4x4 &out)
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

    out = ybi::Float4x4(
        r00, r01, r02, itx, r10, r11, r12, ity, r20, r21, r22, itz, 0.0f, 0.0f, 0.0f, 1.0f);
    return true;
}

static ybi::Float4x4 BuildRasterFromNdc(int width, int height)
{
    const float sx = 0.5f * static_cast<float>(width);
    const float sy = 0.5f * static_cast<float>(height);
    return ybi::Float4x4(
        sx, 0.0f, 0.0f, sx, 0.0f, -sy, 0.0f, sy, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f);
}

static bool ComputeMinimumCameraDifferentials(RenderCameraOverride *camera)
{
    if (!camera || camera->width <= 0 || camera->height <= 0)
    {
        return false;
    }

    LaunchParams params = {};
    params.width = camera->width;
    params.height = camera->height;
    params.cameraOrigin = camera->origin;
    params.cameraFromWorld = camera->cameraFromWorld;
    params.worldFromCamera = camera->worldFromCamera;
    params.cameraFromRaster = camera->cameraFromRaster;
    params.rasterFromCamera = camera->rasterFromCamera;

    bool foundSample = false;
    float minPosXLen = std::numeric_limits<float>::max();
    float minPosYLen = std::numeric_limits<float>::max();
    float minDirXLen = std::numeric_limits<float>::max();
    float minDirYLen = std::numeric_limits<float>::max();
    ybi::Vec3 minPosX = ybi::Vec3(0.0f);
    ybi::Vec3 minPosY = ybi::Vec3(0.0f);
    ybi::Vec3 minDirX = ybi::Vec3(0.0f);
    ybi::Vec3 minDirY = ybi::Vec3(0.0f);
    const int sampleCount = 512;
    const float sampleMaxX = static_cast<float>(std::max(camera->width - 1, 0));
    const float sampleMaxY = static_cast<float>(std::max(camera->height - 1, 0));

    for (int i = 0; i < sampleCount; ++i)
    {
        const float t = sampleCount > 1 ? static_cast<float>(i) / static_cast<float>(sampleCount - 1)
                                        : 0.0f;
        const float sampleX = t * sampleMaxX;
        const float sampleY = t * sampleMaxY;
        const ybi::render::integrator::RayDifferential rayDiff =
            ybi::render::integrator::InitPerspectiveRayDifferential(
                params,
                sampleX,
                sampleY,
                static_cast<unsigned int>(camera->width),
                static_cast<unsigned int>(camera->height));
        if (!rayDiff.valid)
        {
            continue;
        }

        foundSample = true;
        const ybi::Vec3 dox =
            ybi::TransformVectorAffine(camera->cameraFromWorld, rayDiff.originX - rayDiff.origin);
        const ybi::Vec3 doy =
            ybi::TransformVectorAffine(camera->cameraFromWorld, rayDiff.originY - rayDiff.origin);
        const float doxLen = ybi::Length(dox);
        const float doyLen = ybi::Length(doy);
        if (doxLen < minPosXLen)
        {
            minPosXLen = doxLen;
            minPosX = dox;
        }
        if (doyLen < minPosYLen)
        {
            minPosYLen = doyLen;
            minPosY = doy;
        }

        const ybi::Vec3 dir = ybi::Normalize(rayDiff.dir);
        const ybi::Vec3 dirX = ybi::Normalize(rayDiff.dirX);
        const ybi::Vec3 dirY = ybi::Normalize(rayDiff.dirY);
        const ybi::Float4x4 localFromWorld =
            ybi::render::integrator::RotateFromTo(dir, ybi::Vec3(0.0f, 0.0f, 1.0f));
        const ybi::Vec3 df = ybi::Normalize(ybi::TransformVectorAffine(localFromWorld, dir));
        const ybi::Vec3 dxf = ybi::Normalize(ybi::TransformVectorAffine(localFromWorld, dirX));
        const ybi::Vec3 dyf = ybi::Normalize(ybi::TransformVectorAffine(localFromWorld, dirY));
        const ybi::Vec3 ddx = dxf - df;
        const ybi::Vec3 ddy = dyf - df;
        const float ddxLen = ybi::Length(ddx);
        const float ddyLen = ybi::Length(ddy);
        if (ddxLen < minDirXLen)
        {
            minDirXLen = ddxLen;
            minDirX = ddx;
        }
        if (ddyLen < minDirYLen)
        {
            minDirYLen = ddyLen;
            minDirY = ddy;
        }
    }

    if (!foundSample)
    {
        return false;
    }

    camera->minPosDifferentialX = minPosX;
    camera->minPosDifferentialY = minPosY;
    camera->minDirDifferentialX = minDirX;
    camera->minDirDifferentialY = minDirY;
    return ybi::render::integrator::IsFiniteVec3(camera->minPosDifferentialX) &&
           ybi::render::integrator::IsFiniteVec3(camera->minPosDifferentialY) &&
           ybi::render::integrator::IsFiniteVec3(camera->minDirDifferentialX) &&
           ybi::render::integrator::IsFiniteVec3(camera->minDirDifferentialY);
}

static bool FinalizeRenderCameraOverride(RenderCameraOverride *camera)
{
    if (!camera || camera->width <= 0 || camera->height <= 0)
    {
        return false;
    }

    if (!InvertAffine(camera->cameraFromWorld, camera->worldFromCamera))
    {
        return false;
    }
    camera->rasterFromCamera =
        BuildRasterFromNdc(camera->width, camera->height) * camera->clipFromCamera;
    if (!ybi::Invert(camera->rasterFromCamera, &camera->cameraFromRaster))
    {
        return false;
    }

    camera->origin = ybi::TransformPointAffine(camera->worldFromCamera, ybi::Vec3(0.0f));
    return ComputeMinimumCameraDifferentials(camera);
}

static std::optional<RenderCameraOverride>
BuildFallbackRenderCameraOverride(const ybi::Vec3 &eye,
                                  const ybi::Vec3 &lookAt,
                                  int width,
                                  int height,
                                  ybi::UpAxis upAxis)
{
    ybi::Vec3 forward = ybi::Normalize(lookAt - eye);
    if (ybi::Length(forward) <= 1.0e-8f)
    {
        forward = ybi::Vec3(0.0f, 0.0f, 1.0f);
    }
    const ybi::Vec3 worldUp = ybi::ResolveCameraWorldUp(forward, upAxis);
    const ybi::Vec3 right = ybi::Normalize(ybi::Cross(forward, worldUp));
    const ybi::Vec3 up = ybi::Normalize(ybi::Cross(right, forward));
    const float aspect = static_cast<float>(width) / static_cast<float>(std::max(height, 1));
    const float fovY = 45.0f * ybi::kDegToRad;
    const float tanHalfFov = std::tan(fovY * 0.5f);

    RenderCameraOverride camera = {};
    camera.width = width;
    camera.height = height;
    camera.cameraFromWorld = ybi::BuildCameraFromWorld(eye, lookAt, upAxis);
    camera.clipFromCamera = ybi::BuildPerspectiveClipFromCamera(45.0f, width, height);
    camera.U = right * (aspect * tanHalfFov);
    camera.V = up * tanHalfFov;
    camera.W = forward;
    if (!FinalizeRenderCameraOverride(&camera))
    {
        return std::nullopt;
    }
    return camera;
}

static void PopulateLaunchCameraParams(LaunchParams *params, const RenderCameraOverride &camera)
{
    YBI_ASSERT(params);
    params->cameraOrigin = camera.origin;
    params->cameraFromWorld = camera.cameraFromWorld;
    params->worldFromCamera = camera.worldFromCamera;
    params->cameraFromRaster = camera.cameraFromRaster;
    params->rasterFromCamera = camera.rasterFromCamera;
    params->minPosDifferentialX = camera.minPosDifferentialX;
    params->minPosDifferentialY = camera.minPosDifferentialY;
    params->minDirDifferentialX = camera.minDirDifferentialX;
    params->minDirDifferentialY = camera.minDirDifferentialY;
    params->cameraU = camera.U;
    params->cameraV = camera.V;
    params->cameraW = camera.W;
}

static std::optional<RenderCameraOverride> BuildUsdRenderCamera(const Camera &camera)
{
    if (camera.viewportWidth <= 0 || camera.viewportHeight <= 0)
    {
        return std::nullopt;
    }

    ybi::Float4x4 worldFromCamera = {};
    if (!InvertAffine(camera.cameraFromWorld, worldFromCamera))
    {
        return std::nullopt;
    }

    ybi::Vec3 right = ybi::Normalize(ybi::Vec3(
        worldFromCamera.m[0][0], worldFromCamera.m[1][0], worldFromCamera.m[2][0]));
    ybi::Vec3 up = ybi::Normalize(ybi::Vec3(
        worldFromCamera.m[0][1], worldFromCamera.m[1][1], worldFromCamera.m[2][1]));
    ybi::Vec3 forward = ybi::Normalize(ybi::Vec3(
        -worldFromCamera.m[0][2], -worldFromCamera.m[1][2], -worldFromCamera.m[2][2]));
    if (std::abs(ybi::Dot(forward, forward)) < 1e-6f)
    {
        return std::nullopt;
    }

    const float m00 = camera.clipFromCamera.m[0][0];
    const float m11 = camera.clipFromCamera.m[1][1];
    float tanHalfFov = std::tan(45.0f * 0.5f * ybi::kDegToRad);
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
    out.cameraFromWorld = camera.cameraFromWorld;
    out.clipFromCamera = camera.clipFromCamera;
    out.U = right * (aspect * tanHalfFov);
    out.V = up * tanHalfFov;
    out.W = forward;
    if (!FinalizeRenderCameraOverride(&out))
    {
        return std::nullopt;
    }
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
    options.viewportWidth = std::max(1, camera.viewportWidth);
    options.viewportHeight = std::max(1, camera.viewportHeight);
    options.cameraFromWorld = camera.cameraFromWorld;
    options.clipFromCamera = camera.clipFromCamera;

    const size_t srcCount = rootScene->subdivisionMeshes.size();
    const size_t dstBase = rootScene->meshes.size();
    rootScene->meshes.resize(dstBase + srcCount);
    std::atomic<size_t> failureIndex(srcCount);
    tbb::parallel_for(size_t(0), srcCount, [&](size_t i) {
        if (failureIndex.load(std::memory_order_relaxed) != srcCount)
        {
            return;
        }

        const SubdivisionMesh &subdivMesh = rootScene->subdivisionMeshes[i];
        SubdivisionRunResult result = {};
        if (!SubdivideAdaptive(subdivMesh, options, &result))
        {
            size_t expected = srcCount;
            (void)failureIndex.compare_exchange_strong(
                expected, i, std::memory_order_relaxed, std::memory_order_relaxed);
            return;
        }

        result.mesh.materialIndex = subdivMesh.materialIndex;
        result.mesh.parentFromLocal = subdivMesh.parentFromLocal;
        rootScene->meshes[dstBase + i] = std::move(result.mesh);
    });

    const size_t failedIndex = failureIndex.load(std::memory_order_relaxed);
    if (failedIndex != srcCount)
    {
        rootScene->meshes.resize(dstBase);
        fprintf(stderr,
                "Subdivision failed for rootScene->subdivisionMeshes[%zu] prim=%s\n",
                failedIndex,
                rootScene->subdivisionMeshes[failedIndex].primPath.c_str());
        return false;
    }
    rootScene->subdivisionMeshes.clear();
    return true;
}

static bool ParseFloat3(int argc, char **argv, int startIndex, ybi::Vec3 &valueOut)
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
            if (std::find(outPurposes->begin(), outPurposes->end(), normalized) ==
                outPurposes->end())
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

static bool ParseTextureViewSemantic(const std::string &input,
                                     MaterialTextureSemantic *outSemantic)
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
    if (token == "specular" || token == "specularcolor")
    {
        *outSemantic = MaterialTextureSemantic::SpecularColor;
        return true;
    }
    if (token == "clearcoat" || token == "coat")
    {
        *outSemantic = MaterialTextureSemantic::Clearcoat;
        return true;
    }
    if (token == "clearcoatroughness" || token == "coatroughness")
    {
        *outSemantic = MaterialTextureSemantic::ClearcoatRoughness;
        return true;
    }
    return false;
}

static bool TryMapInputNameToSemantic(const std::string &inputName,
                                      MaterialTextureSemantic *outSemantic)
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
    if (inputName == "specularColor")
    {
        *outSemantic = MaterialTextureSemantic::SpecularColor;
        return true;
    }
    if (inputName == "clearcoat")
    {
        *outSemantic = MaterialTextureSemantic::Clearcoat;
        return true;
    }
    if (inputName == "clearcoatRoughness")
    {
        *outSemantic = MaterialTextureSemantic::ClearcoatRoughness;
        return true;
    }
    return false;
}

static void PrintUsage(const char *exeName)
{
    printf("Usage: %s [--file path] [--out path] "
           "[--integrator primary|ao|path] [--spp N] [--max-depth N] "
           "[--device gpu|cpu] [--cam-pos x y z] [--look-at x y z] [--camera path-or-name] "
           "[--ntc] [--view name] "
           "[--purposes csv] [--purpose name]\n",
           exeName);
    printf("  --file USDA/USD path\n");
    printf("  --out PNG output path\n");
    printf("  --integrator primary|ao|path\n");
    printf("  --spp spp passes\n");
    printf("  --max-depth max path depth for --integrator path\n");
    printf("  --device gpu|cpu (cpu requires Embree)\n");
    printf("  --cam-pos optional camera position override\n");
    printf("  --look-at optional look-at target (default bounds center)\n");
    printf("  --camera optional USD camera prim path or unique suffix/name\n");
    printf("  --ntc enable USD NTC decode path (falls back to image textures)\n");
    printf("  --virtual-texture run feedback prepass raygen before beauty\n");
    printf("  --vt-tiles-dir path to *.tiles.bin directory for --virtual-texture\n");
    printf("  --vt-cache-bytes byte budget for stream-page physical cache\n");
    printf("  --vt-tail-max-dim pin mips where max(width,height)<=N into tail cache\n");
    printf("  --vt-max-page-uploads max new stream pages uploaded after each feedback pass (0=unlimited)\n");
    printf("  --pixel x y render only image-space pixel (x right, y down)\n");
    printf("  --write-feedback dump feedback text files next to output PNG\n");
    printf("  --view diffuse|roughness|metallic|occlusion|normal|ior|emissive|opacity|specular|clearcoat|clearcoatroughness\n");
    printf("  --purposes comma-separated default,render,proxy,guide\n");
    printf("  --purpose single purpose token, repeatable; overrides defaults\n");
}

static CliOptions ParseCli(int argc, char **argv)
{
    CliOptions options = {};
    bool purposeOverrideSet = false;
    bool virtualTextureCacheBytesExplicit = false;
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
            else if (value == "path")
            {
                options.integrator = IntegratorType::Path;
            }
            else
            {
                PrintUsage(argv[0]);
                std::abort();
            }
            continue;
        }
        if (arg == "--max-depth")
        {
            if (i + 1 >= argc)
            {
                PrintUsage(argv[0]);
                std::abort();
            }
            options.maxDepth = std::max(0, std::stoi(argv[++i]));
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
        if (arg == "--device")
        {
            if (i + 1 >= argc)
            {
                PrintUsage(argv[0]);
                std::abort();
            }
            const std::string value = argv[++i];
            if (value == "gpu")
            {
                options.deviceKind = DeviceKind::GPU;
            }
            else if (value == "cpu")
            {
                options.deviceKind = DeviceKind::CPU;
            }
            else
            {
                PrintUsage(argv[0]);
                std::abort();
            }
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
            ybi::Vec3 value = {};
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
            ybi::Vec3 value = {};
            if (!ParseFloat3(argc, argv, i + 1, value))
            {
                PrintUsage(argv[0]);
                std::abort();
            }
            options.lookAt = value;
            i += 3;
            continue;
        }
        if (arg == "--camera")
        {
            if (i + 1 >= argc)
            {
                PrintUsage(argv[0]);
                std::abort();
            }
            options.usdCamera = argv[++i];
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
        if (arg == "--virtual-texture")
        {
            options.virtualTexture = true;
            continue;
        }
        if (arg == "--vt-tiles-dir")
        {
            if (i + 1 >= argc)
            {
                PrintUsage(argv[0]);
                std::abort();
            }
            options.virtualTextureTilesDir = argv[++i];
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
        if (arg == "--vt-cache-bytes")
        {
            if (i + 1 >= argc)
            {
                PrintUsage(argv[0]);
                std::abort();
            }
            options.virtualTextureCacheBytes = static_cast<uint64_t>(std::stoull(argv[++i]));
            virtualTextureCacheBytesExplicit = true;
            continue;
        }
        if (arg == "--vt-tail-max-dim")
        {
            if (i + 1 >= argc)
            {
                PrintUsage(argv[0]);
                std::abort();
            }
            options.virtualTextureTailMaxDim = std::max(1, std::stoi(argv[++i]));
            continue;
        }
        if (arg == "--vt-max-page-uploads")
        {
            if (i + 1 >= argc)
            {
                PrintUsage(argv[0]);
                std::abort();
            }
            options.virtualTextureMaxPageUploads = std::max(0, std::stoi(argv[++i]));
            continue;
        }
        if (arg == "--pixel")
        {
            if (i + 2 >= argc)
            {
                PrintUsage(argv[0]);
                std::abort();
            }
            options.singlePixelMode = true;
            options.singlePixelX = std::stoi(argv[++i]);
            options.singlePixelY = std::stoi(argv[++i]);
            continue;
        }
        if (arg == "--write-feedback")
        {
            options.writeFeedbackFiles = true;
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

    if (!virtualTextureCacheBytesExplicit)
    {
        options.virtualTextureCacheBytes =
            options.deviceKind == DeviceKind::CPU ? (8ull << 30u) : (1ull << 30u);
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
    DeviceTextureFormat format = DeviceTextureFormat::RGBA8_UNORM;
    std::vector<uint8_t> pixels;
    std::string sourcePath;
    TextureWrapMode wrapS = TEXTURE_WRAP_MODE_REPEAT;
    TextureWrapMode wrapT = TEXTURE_WRAP_MODE_REPEAT;
};

struct ImageDecodeStats
{
    double totalMs = 0.0;
    double udimResolveMs = 0.0;
    double exrLoadMs = 0.0;
    double exrConvertMs = 0.0;
    double stbiLoadMs = 0.0;
    double stbiCopyMs = 0.0;
    size_t decodedBytes = 0;
    int decodedTiles = 0;
    int cacheHits = 0;
    int cacheMisses = 0;
    int exrTiles = 0;
    int stbiTiles = 0;
    std::vector<std::pair<double, std::string>> slowTiles;
};

static void RecordSlowTile(ImageDecodeStats *stats, double ms, const std::string &path)
{
    if (!stats)
    {
        return;
    }
    if (stats->slowTiles.size() < 10)
    {
        stats->slowTiles.emplace_back(ms, path);
    }
    else
    {
        auto minIt = std::min_element(stats->slowTiles.begin(), stats->slowTiles.end());
        if (minIt != stats->slowTiles.end() && minIt->first < ms)
        {
            *minIt = {ms, path};
        }
    }
    std::sort(stats->slowTiles.begin(), stats->slowTiles.end(),
              [](const auto &a, const auto &b) { return a.first > b.first; });
}

struct UploadedMaterialTextures
{
    std::vector<LaunchParams::MaterialTextureRef> refs;
    std::vector<DeviceTexture> textures;
    uint64_t deviceBytes = 0u;
};

struct UploadedStandaloneTexture
{
    LaunchParams::MaterialTextureRef ref = {};
    DeviceTexture texture = {};
    uint64_t deviceBytes = 0u;
};

static constexpr uint32_t kUdimMin = 1001u;
static constexpr uint32_t kUdimMax = 1100u;
static constexpr int kUdimSlotCount = 128;
static constexpr int kMaterialSemanticCount = static_cast<int>(MaterialTextureSemantic::Count);

static int
MaterialTextureSlotIndex(size_t materialIndex, MaterialTextureSemantic semantic, uint32_t udim)
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

static const MaterialTextureInput *FindTextureInputBySemantic(const MaterialInfo &material,
                                                              MaterialTextureSemantic semantic)
{
    for (const MaterialTextureInput &input : material.textureInputs)
    {
        if (input.texturePath.empty())
        {
            continue;
        }
        MaterialTextureSemantic inputSemantic = MaterialTextureSemantic::Diffuse;
        if (TryMapInputNameToSemantic(input.inputName, &inputSemantic) &&
            inputSemantic == semantic)
        {
            return &input;
        }
    }
    return nullptr;
}

static bool LoadExrRgba8(const std::string &path,
                         int *outWidth,
                         int *outHeight,
                         std::vector<uint8_t> *outRgba8,
                         std::string *outReason,
                         bool flipVertical)
{
    YBI_ASSERT(outWidth);
    YBI_ASSERT(outHeight);
    YBI_ASSERT(outRgba8);

    *outWidth = 0;
    *outHeight = 0;
    outRgba8->clear();

    std::vector<float> rgba;
    if (!ybi::texture::LoadExrRgba(path, outWidth, outHeight, &rgba, outReason, flipVertical))
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

static float ClampFiniteDomeHdrValue(float value)
{
    if (!(value == value) || value <= 0.0f)
    {
        return 0.0f;
    }
    return std::min(value, 65504.0f);
}

static bool LoadExrRgba16Float(const std::string &path,
                               int *outWidth,
                               int *outHeight,
                               std::vector<uint8_t> *outPixels,
                               std::string *outReason,
                               bool flipVertical)
{
    YBI_ASSERT(outWidth);
    YBI_ASSERT(outHeight);
    YBI_ASSERT(outPixels);

    *outWidth = 0;
    *outHeight = 0;
    outPixels->clear();

    std::vector<float> rgba;
    if (!ybi::texture::LoadExrRgba(path, outWidth, outHeight, &rgba, outReason, flipVertical))
    {
        return false;
    }

    const size_t sampleCount = static_cast<size_t>(*outWidth) * static_cast<size_t>(*outHeight) * 4u;
    outPixels->resize(sampleCount * sizeof(uint16_t));
    for (size_t i = 0; i < sampleCount; ++i)
    {
        const uint16_t halfBits = ybi::util::FloatToHalfBits(ClampFiniteDomeHdrValue(rgba[i]));
        std::memcpy(outPixels->data() + i * sizeof(uint16_t), &halfBits, sizeof(halfBits));
    }
    return true;
}

static bool LoadImageRgba8(const std::string &path,
                           int *outWidth,
                           int *outHeight,
                           std::vector<uint8_t> *outRgba8,
                           std::string *outReason,
                           bool flipVertical)
{
    YBI_ASSERT(outWidth);
    YBI_ASSERT(outHeight);
    YBI_ASSERT(outRgba8);

    *outWidth = 0;
    *outHeight = 0;
    outRgba8->clear();

    if (ybi::texture::LowerExt(path) == ".exr")
    {
        return LoadExrRgba8(path, outWidth, outHeight, outRgba8, outReason, flipVertical);
    }

    int channels = 0;
    stbi_set_flip_vertically_on_load(flipVertical ? 1 : 0);
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
    using Clock = std::chrono::steady_clock;
    auto ElapsedMs = [](Clock::time_point start, Clock::time_point end) {
        return std::chrono::duration<double, std::milli>(end - start).count();
    };

    ImageDecodeStats stats = {};
    double cacheHitMs = 0.0;
    const auto totalStart = Clock::now();

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

            std::vector<std::pair<uint32_t, std::string>> udimPaths;
            std::string udimReason;
            const auto udimStart = Clock::now();
            if (!ybi::usd_ntc::CollectUdimPaths(input.texturePath, udimPaths, udimReason))
            {
                std::printf(
                    "Image runtime: failed to resolve UDIMs for material %zu input %s (%s): %s\n",
                    materialIndex,
                    input.inputName.c_str(),
                    input.texturePath.c_str(),
                    udimReason.c_str());
                stats.udimResolveMs += ElapsedMs(udimStart, Clock::now());
                continue;
            }
            stats.udimResolveMs += ElapsedMs(udimStart, Clock::now());

            for (const auto &entry : udimPaths)
            {
                const uint32_t udim = entry.first;
                const std::string &tilePath = entry.second;
                const int dstIndex = MaterialTextureSlotIndex(materialIndex, semantic, udim);
                YBI_ASSERT(dstIndex >= 0 && static_cast<size_t>(dstIndex) < outTextures->size());

                auto cached = decodedByPath.find(tilePath);
                if (cached != decodedByPath.end())
                {
                    const auto cacheStart = Clock::now();
                    DecodedMaterialTexture tile = cached->second;
                    tile.udim = udim;
                    tile.wrapS = input.wrapS;
                    tile.wrapT = input.wrapT;
                    (*outTextures)[static_cast<size_t>(dstIndex)] = std::move(tile);
                    decodedTiles++;
                    stats.cacheHits++;
                    stats.decodedTiles++;
                    cacheHitMs += ElapsedMs(cacheStart, Clock::now());
                    continue;
                }
                stats.cacheMisses++;

                DecodedMaterialTexture texture = {};
                texture.valid = true;
                texture.udim = udim;
                texture.sourcePath = tilePath;
                texture.wrapS = input.wrapS;
                texture.wrapT = input.wrapT;

                std::string reason;
                const auto loadStart = Clock::now();
                if (!LoadImageRgba8(tilePath,
                                    &texture.width,
                                    &texture.height,
                                    &texture.pixels,
                                    &reason,
                                    true))
                {
                    std::printf(
                        "Image runtime: failed to load material %zu input %s tile %u (%s): %s\n",
                        materialIndex,
                        input.inputName.c_str(),
                        udim,
                        tilePath.c_str(),
                        reason.c_str());
                    continue;
                }
                const double loadMs = ElapsedMs(loadStart, Clock::now());
                const bool isExr = (ybi::texture::LowerExt(tilePath) == ".exr");
                if (isExr)
                {
                    stats.exrLoadMs += loadMs;
                    stats.exrTiles++;
                }
                else
                {
                    stats.stbiLoadMs += loadMs;
                    stats.stbiTiles++;
                }
                stats.decodedBytes += texture.pixels.size();
                stats.decodedTiles++;
                RecordSlowTile(&stats, loadMs, tilePath);

                (*outTextures)[static_cast<size_t>(dstIndex)] = texture;
                decodedByPath.emplace(tilePath, std::move(texture));
                decodedTiles++;
            }
        }
    }

    stats.totalMs = ElapsedMs(totalStart, Clock::now());

    std::printf("Image runtime: decoded %d/%zu material UDIM tiles\n",
                decodedTiles,
                materials.size() * static_cast<size_t>(kMaterialSemanticCount) *
                    static_cast<size_t>(kUdimSlotCount));
    const double imageLoadMs = stats.exrLoadMs + stats.stbiLoadMs;
    const double otherMs =
        std::max(0.0, stats.totalMs - (stats.udimResolveMs + imageLoadMs + cacheHitMs));
    const char *bottleneck = "other";
    double bottleneckMs = otherMs;
    if (stats.udimResolveMs > bottleneckMs)
    {
        bottleneck = "udim_resolve";
        bottleneckMs = stats.udimResolveMs;
    }
    if (cacheHitMs > bottleneckMs)
    {
        bottleneck = "cache_hit_copy";
        bottleneckMs = cacheHitMs;
    }
    if (imageLoadMs > bottleneckMs)
    {
        bottleneck = "image_load";
        bottleneckMs = imageLoadMs;
    }
    std::printf(
        "Image runtime: decode timing total %.2f ms (udim %.2f ms, cache-hit %.2f ms, "
        "image-load %.2f ms [exr %.2f ms/%d, stbi %.2f ms/%d], other %.2f ms)\n",
        stats.totalMs,
        stats.udimResolveMs,
        cacheHitMs,
        imageLoadMs,
        stats.exrLoadMs,
        stats.exrTiles,
        stats.stbiLoadMs,
        stats.stbiTiles,
        otherMs);
    std::printf("Image runtime: bottleneck guess %s (%.2f ms)\n", bottleneck, bottleneckMs);
    std::printf("Image runtime: cache hits %d misses %d decoded bytes %zu\n",
                stats.cacheHits,
                stats.cacheMisses,
                stats.decodedBytes);
    if (!stats.slowTiles.empty())
    {
        std::printf("Image runtime: slow tiles (ms):\n");
        for (const auto &entry : stats.slowTiles)
        {
            std::printf("  %.2f %s\n", entry.first, entry.second.c_str());
        }
    }
    return true;
}

static DeviceTextureWrapMode ToDeviceTextureWrapMode(const ybi::TextureWrapMode mode)
{
    switch (mode)
    {
        case TEXTURE_WRAP_MODE_REPEAT:
            return DeviceTextureWrapMode::Repeat;
        case TEXTURE_WRAP_MODE_CLAMP:
            return DeviceTextureWrapMode::Clamp;
        case TEXTURE_WRAP_MODE_MIRROR:
            return DeviceTextureWrapMode::Mirror;
        case TEXTURE_WRAP_MODE_BLACK:
            return DeviceTextureWrapMode::Black;
        case TEXTURE_WRAP_MODE_USE_METADATA:
            return DeviceTextureWrapMode::UseMetadata;
        case TEXTURE_WRAP_MODE_UNKNOWN:
        default:
            return DeviceTextureWrapMode::Unknown;
    }
}

static bool UploadDecodedTextures(Device *device,
                                  const std::vector<DecodedMaterialTexture> &decodedTextures,
                                  UploadedMaterialTextures *outTextures,
                                  std::string *outError)
{
    YBI_ASSERT(device);
    YBI_ASSERT(outTextures);
    for (DeviceTexture &texture : outTextures->textures)
    {
        device->DestroyTexture(texture);
    }
    outTextures->textures.clear();
    outTextures->refs.clear();
    outTextures->deviceBytes = 0u;
    outTextures->refs.resize(decodedTextures.size());

    for (size_t i = 0; i < decodedTextures.size(); ++i)
    {
        const DecodedMaterialTexture &src = decodedTextures[i];
        if (!src.valid || src.width <= 0 || src.height <= 0 || src.pixels.empty())
        {
            continue;
        }

        DeviceTextureCreateInfo createInfo = {};
        createInfo.pixels = src.pixels.data();
        createInfo.pixelBytes = src.pixels.size();
        createInfo.width = static_cast<uint32_t>(src.width);
        createInfo.height = static_cast<uint32_t>(src.height);
        createInfo.wrapS = ToDeviceTextureWrapMode(src.wrapS);
        createInfo.wrapT = ToDeviceTextureWrapMode(src.wrapT);
        createInfo.filter = DeviceTextureFilterMode::Nearest;
        createInfo.format = src.format;

        DeviceTexture texture = {};
        if (!device->CreateTexture(createInfo, &texture, outError))
        {
            return false;
        }

        outTextures->textures.push_back(texture);
        outTextures->deviceBytes += static_cast<uint64_t>(src.pixels.size());
        outTextures->refs[i].textureObject = texture.handle;
        outTextures->refs[i].width = src.width;
        outTextures->refs[i].height = src.height;
        outTextures->refs[i].valid = 1;
        outTextures->refs[i].wrapS = static_cast<int>(src.wrapS);
        outTextures->refs[i].wrapT = static_cast<int>(src.wrapT);
        outTextures->refs[i].format = src.format;
    }

    return true;
}

static void DestroyUploadedTextures(Device *device, UploadedMaterialTextures *textures)
{
    YBI_ASSERT(device);
    if (!textures)
    {
        return;
    }
    for (DeviceTexture &texture : textures->textures)
    {
        device->DestroyTexture(texture);
    }
    textures->textures.clear();
    textures->refs.clear();
    textures->deviceBytes = 0u;
}

static bool UploadDecodedTexture(Device *device,
                                 const DecodedMaterialTexture &decodedTexture,
                                 UploadedStandaloneTexture *outTexture,
                                 std::string *outError)
{
    YBI_ASSERT(device);
    YBI_ASSERT(outTexture);
    if (outTexture->texture.valid)
    {
        device->DestroyTexture(outTexture->texture);
    }
    outTexture->ref = {};
    outTexture->texture = {};
    outTexture->deviceBytes = 0u;

    if (!decodedTexture.valid || decodedTexture.width <= 0 || decodedTexture.height <= 0 ||
        decodedTexture.pixels.empty())
    {
        return true;
    }

    DeviceTextureCreateInfo createInfo = {};
    createInfo.pixels = decodedTexture.pixels.data();
    createInfo.pixelBytes = decodedTexture.pixels.size();
    createInfo.width = static_cast<uint32_t>(decodedTexture.width);
    createInfo.height = static_cast<uint32_t>(decodedTexture.height);
    createInfo.wrapS = ToDeviceTextureWrapMode(decodedTexture.wrapS);
    createInfo.wrapT = ToDeviceTextureWrapMode(decodedTexture.wrapT);
    createInfo.filter = DeviceTextureFilterMode::Linear;
    createInfo.format = decodedTexture.format;

    if (!device->CreateTexture(createInfo, &outTexture->texture, outError))
    {
        return false;
    }

    outTexture->deviceBytes = static_cast<uint64_t>(decodedTexture.pixels.size());
    outTexture->ref.textureObject = outTexture->texture.handle;
    outTexture->ref.width = decodedTexture.width;
    outTexture->ref.height = decodedTexture.height;
    outTexture->ref.valid = 1;
    outTexture->ref.wrapS = static_cast<int>(decodedTexture.wrapS);
    outTexture->ref.wrapT = static_cast<int>(decodedTexture.wrapT);
    outTexture->ref.format = decodedTexture.format;
    return true;
}

static void DestroyUploadedTexture(Device *device, UploadedStandaloneTexture *texture)
{
    YBI_ASSERT(device);
    if (!texture)
    {
        return;
    }
    if (texture->texture.valid)
    {
        device->DestroyTexture(texture->texture);
    }
    texture->ref = {};
    texture->texture = {};
    texture->deviceBytes = 0u;
}

struct UploadedMeshRefs
{
    std::vector<LaunchParams::InstanceGeomRef> refs;
    std::vector<DeviceMemoryView<uint8_t>> ownedBuffers;
    uint64_t deviceBytes = 0u;
};

static bool PathMatchesOrContainsChild(const std::string &primPath, const std::string &targetPath)
{
    if (primPath.empty() || targetPath.empty())
    {
        return false;
    }
    if (primPath == targetPath)
    {
        return true;
    }
    return primPath.size() > targetPath.size() &&
           primPath.compare(0, targetPath.size(), targetPath) == 0 &&
           primPath[targetPath.size()] == '/';
}

static void ResolveLightShadowExcludes(const std::vector<SceneMeshUploadRef> &meshUploadRefs,
                                       const std::vector<LightInfo> &lights,
                                       std::vector<PackedLight> *outPackedLights,
                                       std::vector<int> *outExcludeRefs)
{
    YBI_ASSERT(outPackedLights);
    YBI_ASSERT(outExcludeRefs);
    outPackedLights->clear();
    outExcludeRefs->clear();
    outPackedLights->reserve(lights.size());

    for (const LightInfo &light : lights)
    {
        PackedLight packed = light.packed;
        packed.shadowExcludeOffset = static_cast<uint32_t>(outExcludeRefs->size());
        packed.shadowExcludeCount = 0u;

        std::unordered_set<uint32_t> seenRefs;
        for (const std::string &targetPath : light.shadowExcludePaths)
        {
            for (const SceneMeshUploadRef &job : meshUploadRefs)
            {
                if (!job.mesh || job.mesh->primPath.empty() ||
                    !PathMatchesOrContainsChild(job.mesh->primPath, targetPath))
                {
                    continue;
                }
                if (seenRefs.insert(job.refIndex).second)
                {
                    outExcludeRefs->push_back(static_cast<int>(job.refIndex));
                }
            }
        }

        packed.shadowExcludeCount =
            static_cast<uint32_t>(outExcludeRefs->size()) - packed.shadowExcludeOffset;
        outPackedLights->push_back(packed);
    }
}

static bool FindFirstDomeTexturePath(const std::vector<LightInfo> &lights, std::string *outPath)
{
    YBI_ASSERT(outPath);
    outPath->clear();
    for (const LightInfo &light : lights)
    {
        if (light.packed.type != static_cast<uint32_t>(LightType::Dome) || light.texturePath.empty())
        {
            continue;
        }
        *outPath = light.texturePath;
        return true;
    }
    return false;
}

struct HarnessMemoryStats
{
    uint64_t hostSourceSceneBytes = 0u;
    uint64_t hostFlattenedSceneBytes = 0u;
    uint64_t hostTessellatedMeshBytes = 0u;
    uint64_t hostDecodedTextureBytes = 0u;
    uint64_t hostMaterialTextureRefBytes = 0u;
    uint64_t virtualTextureVirtualBytes = 0u;
    uint64_t deviceMaterialBytes = 0u;
    uint64_t deviceLightBytes = 0u;
    uint64_t deviceMaterialTextureRefBytes = 0u;
    uint64_t deviceUploadedTextureBytes = 0u;
    uint64_t deviceMeshUploadBytes = 0u;
    uint64_t deviceInstanceRefBytes = 0u;
    uint64_t deviceBvhBytes = 0u;
};

static double BytesToMiB(uint64_t bytes)
{
    return double(bytes) / (1024.0 * 1024.0);
}

static void PrintByteStat(const char *label, uint64_t bytes)
{
    std::printf("memory: %-28s %12llu bytes  %10.2f MiB\n",
                label,
                static_cast<unsigned long long>(bytes),
                BytesToMiB(bytes));
}

static void PrintSetupMemoryStats(const HarnessMemoryStats &stats)
{
    const uint64_t totalHost = stats.hostSourceSceneBytes + stats.hostFlattenedSceneBytes +
                               stats.hostDecodedTextureBytes + stats.hostMaterialTextureRefBytes;
    const uint64_t totalDevice = stats.deviceMaterialBytes + stats.deviceLightBytes +
                                 stats.deviceMaterialTextureRefBytes + stats.deviceUploadedTextureBytes +
                                 stats.deviceMeshUploadBytes + stats.deviceInstanceRefBytes +
                                 stats.deviceBvhBytes;
    std::printf("memory: setup\n");
    PrintByteStat("host.scene.source", stats.hostSourceSceneBytes);
    PrintByteStat("host.scene.flattened", stats.hostFlattenedSceneBytes);
    PrintByteStat("host.scene.tessellated", stats.hostTessellatedMeshBytes);
    PrintByteStat("host.textures.decoded", stats.hostDecodedTextureBytes);
    PrintByteStat("host.texture_refs", stats.hostMaterialTextureRefBytes);
    if (stats.virtualTextureVirtualBytes != 0u)
    {
        PrintByteStat("vt.virtual_source", stats.virtualTextureVirtualBytes);
    }
    PrintByteStat("device.materials", stats.deviceMaterialBytes);
    PrintByteStat("device.lights", stats.deviceLightBytes);
    PrintByteStat("device.texture_refs", stats.deviceMaterialTextureRefBytes);
    PrintByteStat("device.textures.classic", stats.deviceUploadedTextureBytes);
    PrintByteStat("device.mesh_uploads", stats.deviceMeshUploadBytes);
    PrintByteStat("device.instance_refs", stats.deviceInstanceRefBytes);
    PrintByteStat("device.bvh", stats.deviceBvhBytes);
    PrintByteStat("host.total", totalHost);
    PrintByteStat("device.total", totalDevice);
}

static void PrintVirtualTextureMemoryStats(const HarnessMemoryStats &setupStats,
                                           const ybi::texture::VirtualTextureMemoryStats &vtStats,
                                           uint64_t renderDeviceBytes,
                                           uint64_t renderHostBytes)
{
    const uint64_t totalHost = setupStats.hostSourceSceneBytes + setupStats.hostFlattenedSceneBytes +
                               setupStats.hostMaterialTextureRefBytes + vtStats.hostPageTableBytes +
                               vtStats.hostMetaBytes + vtStats.hostTailBytes + vtStats.hostStreamBytes +
                               renderHostBytes;
    const uint64_t totalDevice = setupStats.deviceMaterialTextureRefBytes +
                                 setupStats.deviceMaterialBytes + setupStats.deviceLightBytes +
                                 setupStats.deviceMeshUploadBytes + setupStats.deviceInstanceRefBytes +
                                 setupStats.deviceBvhBytes + vtStats.devicePageTableBytes +
                                 vtStats.deviceMetaBytes + vtStats.deviceTailBytes +
                                 vtStats.deviceStreamBytes + renderDeviceBytes;
    std::printf("memory: virtual-texture resident\n");
    PrintByteStat("vt.host.page_table", vtStats.hostPageTableBytes);
    PrintByteStat("vt.host.meta", vtStats.hostMetaBytes);
    PrintByteStat("vt.host.tail", vtStats.hostTailBytes);
    PrintByteStat("vt.host.stream", vtStats.hostStreamBytes);
    PrintByteStat("vt.device.page_table", vtStats.devicePageTableBytes);
    PrintByteStat("vt.device.meta", vtStats.deviceMetaBytes);
    PrintByteStat("vt.device.tail", vtStats.deviceTailBytes);
    PrintByteStat("vt.device.stream", vtStats.deviceStreamBytes);
    PrintByteStat("render.host", renderHostBytes);
    PrintByteStat("render.device", renderDeviceBytes);
    PrintByteStat("host.total", totalHost);
    PrintByteStat("device.total", totalDevice);
}

static const Attribute *FindMeshSTAttribute(const Mesh &mesh)
{
    for (const Attribute &attr : mesh.attributes)
    {
        if (attr.name == "st" && attr.type == AttributeType::Float2 &&
            attr.interpolation == PrimvarInterpolation::FaceVarying && attr.indices.size() > 0)
        {
            if ((attr.data.size() % sizeof(ybi::Vec2)) == 0)
            {
                return &attr;
            }
        }
    }
    return nullptr;
}

static const Attribute *FindMeshNormalAttribute(const Mesh &mesh)
{
    for (const Attribute &attr : mesh.attributes)
    {
        if (attr.name == "normal" && attr.type == AttributeType::Float3 &&
            attr.interpolation == PrimvarInterpolation::FaceVarying && attr.indices.size() > 0)
        {
            if ((attr.data.size() % sizeof(ybi::Vec3)) == 0)
            {
                return &attr;
            }
        }
    }
    return nullptr;
}

static UploadedMeshRefs
UploadScenePoolMeshRefs(Device *device, const std::vector<SceneMeshUploadRef> &meshUploadRefs)
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

        const size_t positionsBytes = sizeof(ybi::Vec3) * mesh.positions.size();
        const size_t indicesBytes = sizeof(int) * mesh.indices.size();
        DeviceMemoryView<uint8_t> positionsBuffer = device->AllocBytes(positionsBytes);
        DeviceMemoryView<uint8_t> indicesBuffer = device->AllocBytes(indicesBytes);
        out.deviceBytes += positionsBuffer.numBytes() + indicesBuffer.numBytes();
        device->CopyBytesToDevice(positionsBuffer, mesh.positions.data(), positionsBytes);
        device->CopyBytesToDevice(indicesBuffer, mesh.indices.data(), indicesBytes);

        DeviceMemoryView<uint8_t> texcoordsBuffer = {};
        DeviceMemoryView<uint8_t> texcoordIndicesBuffer = {};
        DeviceMemoryView<uint8_t> normalsBuffer = {};
        DeviceMemoryView<uint8_t> normalIndicesBuffer = {};
        int texcoordCount = 0;
        int texcoordIndexCount = 0;
        int normalCount = 0;
        int normalIndexCount = 0;
        const Attribute *stAttr = FindMeshSTAttribute(mesh);
        if (stAttr)
        {
            const size_t texcoordsBytes = stAttr->data.size();
            const size_t texcoordIndicesBytes = sizeof(int) * stAttr->indices.size();
            texcoordsBuffer = device->AllocBytes(texcoordsBytes);
            texcoordIndicesBuffer = device->AllocBytes(texcoordIndicesBytes);
            out.deviceBytes += texcoordsBuffer.numBytes() + texcoordIndicesBuffer.numBytes();
            device->CopyBytesToDevice(texcoordsBuffer, stAttr->data.data(), texcoordsBytes);
            device->CopyBytesToDevice(
                texcoordIndicesBuffer, stAttr->indices.data(), texcoordIndicesBytes);
            out.ownedBuffers.push_back(texcoordsBuffer);
            out.ownedBuffers.push_back(texcoordIndicesBuffer);
            texcoordCount = int(stAttr->data.size() / sizeof(ybi::Vec2));
            texcoordIndexCount = int(stAttr->indices.size());
        }
        const Attribute *normalAttr = FindMeshNormalAttribute(mesh);
        if (normalAttr)
        {
            const size_t normalsBytes = normalAttr->data.size();
            const size_t normalIndicesBytes = sizeof(int) * normalAttr->indices.size();
            normalsBuffer = device->AllocBytes(normalsBytes);
            normalIndicesBuffer = device->AllocBytes(normalIndicesBytes);
            out.deviceBytes += normalsBuffer.numBytes() + normalIndicesBuffer.numBytes();
            device->CopyBytesToDevice(normalsBuffer, normalAttr->data.data(), normalsBytes);
            device->CopyBytesToDevice(
                normalIndicesBuffer, normalAttr->indices.data(), normalIndicesBytes);
            out.ownedBuffers.push_back(normalsBuffer);
            out.ownedBuffers.push_back(normalIndicesBuffer);
            normalCount = int(normalAttr->data.size() / sizeof(ybi::Vec3));
            normalIndexCount = int(normalAttr->indices.size());
        }

        out.refs[job.refIndex] = {(unsigned long long)positionsBuffer.data(),
                                  (unsigned long long)indicesBuffer.data(),
                                  (unsigned long long)texcoordsBuffer.data(),
                                  (unsigned long long)texcoordIndicesBuffer.data(),
                                  (unsigned long long)normalsBuffer.data(),
                                  (unsigned long long)normalIndicesBuffer.data(),
                                  (int)mesh.positions.size(),
                                  (int)mesh.indices.size(),
                                  texcoordCount,
                                  texcoordIndexCount,
                                  normalCount,
                                  normalIndexCount,
                                  mesh.materialIndex};
        out.ownedBuffers.push_back(positionsBuffer);
        out.ownedBuffers.push_back(indicesBuffer);
    }
    return out;
}

static void
RenderTraversable(Device *device,
                  BVHHandle traversable,
                  const ybi::Vec3 &boundsMin,
                  const ybi::Vec3 &boundsMax,
                  const char *outputFile,
                  IntegratorType integrator,
                  int spp,
                  int maxDepth,
                  DevicePtr instanceGeomRefsBuffer,
                  int instanceGeomRefCount,
                  DevicePtr materialsBuffer,
                  int materialCount,
                  DevicePtr lightsBuffer,
                  int lightCount,
                  DevicePtr lightShadowExcludeRefsBuffer,
                  int lightShadowExcludeRefCount,
                  DevicePtr materialTextureRefsBuffer,
                  int materialTextureRefCount,
                  int materialTextureRefStride,
                  int materialTextureRefSemanticCount,
                  const LaunchParams::MaterialTextureRef &domeTextureRef,
                  MaterialTextureSemantic textureViewSemantic,
                  bool virtualTexture,
                  const std::vector<ybi::texture::VirtualTextureRegisterInput> &virtualTextureRegistrations,
                  const HarnessMemoryStats &setupMemoryStats,
                  uint64_t virtualTextureCacheBytes,
                  int virtualTextureTailMaxDim,
                  int virtualTextureMaxPageUploads,
                  bool singlePixelMode,
                  int singlePixelX,
                  int singlePixelY,
                  bool writeFeedbackFiles,
                  ybi::UpAxis fallbackCameraUpAxis,
                  const std::optional<RenderCameraOverride> &cameraOverride,
                  const std::optional<ybi::Vec3> &cameraPositionOverride,
                  const std::optional<ybi::Vec3> &lookAtOverride)
{
    YBI_ASSERT(device);
    using Clock = std::chrono::steady_clock;
    auto LogDuration = [](const char *label, int sppIndex, double ms) {
        std::printf("profile: %s spp=%d %.3f ms\n", label, sppIndex, ms);
    };
    auto LogProcessFeedbackProfile = [](const char *label,
                                        int sppIndex,
                                        const ybi::texture::VirtualTextureUpdateStats &stats) {
        std::printf(
            "profile: %s spp=%d total=%.3f ms histogram=%.3f resolve=%.3f touch=%.3f "
            "allocate=%.3f(%u) load=%.3f(%u) upload=%.3f(%u) update_page_table=%.3f(%u)\n",
            label,
            sppIndex,
            stats.totalMs,
            stats.histogramMs,
            stats.resolveKeyMs,
            stats.touchResidentMs,
            stats.allocateStreamSlotMs,
            stats.allocateStreamSlotCalls,
            stats.loadStreamPageMs,
            stats.loadStreamPageCalls,
            stats.uploadStreamPageMs,
            stats.uploadStreamPageCalls,
            stats.updatePageTableMs,
            stats.updatePageTableCalls);
    };
    printf("render: begin\n");
    fflush(stdout);
    const int width = cameraOverride.has_value() ? cameraOverride->width : 1280;
    const int height = cameraOverride.has_value() ? cameraOverride->height : 720;
    if (singlePixelMode)
    {
        if (singlePixelX < 0 || singlePixelY < 0 || singlePixelX >= width || singlePixelY >= height)
        {
            fprintf(stderr,
                    "single-pixel render: pixel (%d,%d) out of bounds for %dx%d\n",
                    singlePixelX,
                    singlePixelY,
                    width,
                    height);
            std::abort();
        }
        std::printf("single-pixel render: (%d,%d) on %dx%d\n",
                    singlePixelX,
                    singlePixelY,
                    width,
                    height);
    }
    const size_t imageSize = static_cast<size_t>(width) * static_cast<size_t>(height) * 4;
    DeviceMemoryView<uint8_t> imageBuffer = device->AllocBytes(imageSize);
    printf("render: image buffer allocated\n");
    fflush(stdout);

    const int sppPassCount = std::max(1, spp);
    const int feedbackCapacity = std::max(1, width * height);
    const size_t feedbackKeysBytes =
        static_cast<size_t>(feedbackCapacity) * sizeof(unsigned long long);
    const size_t feedbackStatsBytes = 2u * sizeof(unsigned int);
    DeviceMemoryView<uint8_t> feedbackKeysBuffer = device->AllocBytes(feedbackKeysBytes);
    DeviceMemoryView<uint8_t> feedbackStatsBuffer = device->AllocBytes(feedbackStatsBytes);

    std::filesystem::path feedbackDir;
    if (writeFeedbackFiles)
    {
        feedbackDir = std::filesystem::path(outputFile);
        feedbackDir += ".feedback";
        std::error_code feedbackEc;
        std::filesystem::create_directories(feedbackDir, feedbackEc);
        if (feedbackEc)
        {
            fprintf(stderr, "Failed to create feedback dir: %s\n", feedbackDir.string().c_str());
            std::abort();
        }
    }

    const ybi::Vec3 center = (boundsMin + boundsMax) * 0.5f;
    const ybi::Vec3 extent = boundsMax - boundsMin;
    const float diagonal = std::max(0.001f, ybi::Length(extent));

    LaunchParams params = {};
    params.traversable = traversable;
    params.image = reinterpret_cast<DevicePtr>(imageBuffer.data());
    params.width = width;
    params.height = height;
    if (cameraOverride.has_value())
    {
        PopulateLaunchCameraParams(&params, *cameraOverride);
    }
    else
    {
        const ybi::UpAxis upAxis = fallbackCameraUpAxis;
        const ybi::Vec3 eye = cameraPositionOverride.has_value()
                                    ? cameraPositionOverride.value()
                                    : center + ybi::UpAxisVector(upAxis) * (1.25f * diagonal);
        const ybi::Vec3 lookAt = lookAtOverride.has_value() ? lookAtOverride.value() : center;
        const std::optional<RenderCameraOverride> fallbackCamera =
            BuildFallbackRenderCameraOverride(eye, lookAt, width, height, upAxis);
        if (!fallbackCamera.has_value())
        {
            std::fprintf(stderr, "Failed to build fallback camera.\n");
            std::abort();
        }
        PopulateLaunchCameraParams(&params, *fallbackCamera);
    }
    params.wireframe.lineWidth = 0.012f;
    params.wireframe.lineFeather = 0.006f;
    params.wireframe.edgeDarkness = 0.10f;
    params.wireframe.padding = 0.0f;
    params.integrator = integrator == IntegratorType::AO
                            ? 1
                            : (integrator == IntegratorType::Path ? 3 : 0);
    params.spp = std::max(1, spp);
    params.maxDepth = std::max(0, maxDepth);
    params.aoBias = 0.002f * diagonal;
    params.aoMaxDistance = 0.25f * diagonal;
    params.instanceGeomRefs = (unsigned long long)instanceGeomRefsBuffer;
    params.instanceGeomRefCount = instanceGeomRefCount;
    params.materials = (unsigned long long)materialsBuffer;
    params.materialCount = materialCount;
    params.lights = (unsigned long long)lightsBuffer;
    params.lightCount = lightCount;
    params.lightShadowExcludeRefs = (unsigned long long)lightShadowExcludeRefsBuffer;
    params.lightShadowExcludeRefCount = lightShadowExcludeRefCount;
    params.materialTextureRefs = (unsigned long long)materialTextureRefsBuffer;
    params.materialTextureRefCount = materialTextureRefCount;
    params.materialTextureRefStride = materialTextureRefStride;
    params.materialTextureRefSemanticCount = materialTextureRefSemanticCount;
    params.domeTextureRef = domeTextureRef;
    params.textureViewSemantic = static_cast<int>(textureViewSemantic);
    params.feedbackKeys = (unsigned long long)feedbackKeysBuffer.data();
    params.feedbackStats = (unsigned long long)feedbackStatsBuffer.data();
    params.feedbackCapacity = feedbackCapacity;
    params.feedbackSamplePercent = 10;
    params.feedbackTileSize = 128;
    params.currentSpp = 0;
    params.singlePixelEnabled = singlePixelMode ? 1 : 0;
    params.singlePixelX = singlePixelX;
    params.singlePixelY = singlePixelY;
    params.virtualTextureEnabled = virtualTexture ? 1 : 0;
    params.virtualTexturePageTableEntries = 0ull;
    params.virtualTexturePageTableMipOffsets = 0ull;
    params.virtualTexturePageTableMipWidths = 0ull;
    params.virtualTexturePageTableMipHeights = 0ull;
    params.virtualTexturePageTableMipCount = 0;
    params.virtualTexturePageSize = 128;
    params.virtualTextureStreamPixels = 0ull;
    params.virtualTextureStreamPageCountX = 0;
    params.virtualTextureStreamPageCountY = 0;
    params.virtualTextureSampleMip = 0;
    params.virtualTextureTextureMeta = 0ull;
    params.virtualTextureTextureMetaCount = 0;
    params.virtualTextureMipInfos = 0ull;
    params.virtualTextureMipInfoCount = 0;
    params.virtualTextureUdimInfos = 0ull;
    params.virtualTextureUdimInfoCount = 0;

    DeviceMemoryView<uint8_t> paramsBuffer = device->AllocBytes(sizeof(LaunchParams));
    printf("render: params buffer allocated\n");
    fflush(stdout);

    std::vector<uint8_t> hostPassImage(imageSize, 0);
    std::vector<float> accumRgb(static_cast<size_t>(width) * static_cast<size_t>(height) * 3u,
                                0.0f);
    std::vector<unsigned long long> feedbackKeysHost(feedbackCapacity, 0ull);
    unsigned int feedbackStatsHost[2] = {0u, 0u};
    unsigned int feedbackStatsZero[2] = {0u, 0u};
    const uint64_t renderDeviceBytes = imageBuffer.numBytes() + feedbackKeysBuffer.numBytes() +
                                       feedbackStatsBuffer.numBytes() + paramsBuffer.numBytes();
    const uint64_t renderHostBytes = static_cast<uint64_t>(hostPassImage.size()) +
                                     static_cast<uint64_t>(accumRgb.size()) * sizeof(float) +
                                     static_cast<uint64_t>(feedbackKeysHost.size()) *
                                         sizeof(unsigned long long) +
                                     sizeof(feedbackStatsHost) + sizeof(feedbackStatsZero);
    ybi::texture::VirtualTextureManager virtualTextureManager;
    auto WriteFeedbackFile = [&](const char *name,
                                 int sppIndex,
                                 unsigned int sampledCount,
                                 unsigned int copyCount,
                                 unsigned int overflowCount) {
        std::vector<ybi::texture::VirtualTextureFeedbackEntry> histogram;
        if (!writeFeedbackFiles)
        {
            return histogram;
        }
        const std::filesystem::path feedbackPath = feedbackDir / name;
        std::string feedbackError;
        if (!ybi::texture::WriteFeedbackFile(feedbackPath,
                                             sppIndex,
                                             sampledCount,
                                             copyCount,
                                             overflowCount,
                                             feedbackKeysHost,
                                             &histogram,
                                             &feedbackError))
        {
            fprintf(stderr,
                    "Failed to write feedback file: %s (%s)\n",
                    feedbackPath.string().c_str(),
                    feedbackError.c_str());
            std::abort();
        }
        return histogram;
    };

    if (virtualTexture)
    {
        ybi::texture::VirtualTextureManagerConfig vtConfig = {};
        vtConfig.pageSize = 128u;
        vtConfig.cacheBytes = virtualTextureCacheBytes;
        vtConfig.tailMaxDim = static_cast<uint32_t>(std::max(1, virtualTextureTailMaxDim));
        vtConfig.maxUploadsPerPass = static_cast<uint32_t>(std::max(0, virtualTextureMaxPageUploads));
        std::string vtError;
        if (!virtualTextureManager.Initialize(device, vtConfig, &vtError))
        {
            fprintf(stderr, "virtual-texture manager init failed: %s\n", vtError.c_str());
            std::abort();
        }

        for (const ybi::texture::VirtualTextureRegisterInput &reg : virtualTextureRegistrations)
        {
            if (!virtualTextureManager.RegisterTexture(reg, &vtError))
            {
                fprintf(stderr,
                        "virtual-texture manager register failed (material %u): %s\n",
                        reg.textureId,
                        vtError.c_str());
                std::abort();
            }
        }
        if (!virtualTextureManager.Finalize(&vtError))
        {
            fprintf(stderr, "virtual-texture manager finalize failed: %s\n", vtError.c_str());
            std::abort();
        }
        virtualTextureManager.BindLaunchParams(&params);

        params.integrator = 0;
        params.spp = 1;
        params.currentSpp = 0;
        params.feedbackSamplePercent = 100;
        device->CopyBytesToDevice(paramsBuffer, &params, sizeof(LaunchParams));
        device->CopyBytesToDevice(feedbackStatsBuffer, feedbackStatsZero, feedbackStatsBytes);

        DispatchParams dispatchParams = {};
        dispatchParams.width = static_cast<uint32_t>(width);
        dispatchParams.height = static_cast<uint32_t>(height);
        dispatchParams.spp = static_cast<uint32_t>(params.spp);
        dispatchParams.launchParamsDevice = reinterpret_cast<uint64_t>(paramsBuffer.data());
        dispatchParams.launchParamsSize = sizeof(LaunchParams);
        dispatchParams.outputRGBA8 = imageBuffer;
        if (!device->DispatchKernel(RenderKernelId::PrimaryDiffuse, dispatchParams))
        {
            fprintf(stderr, "DispatchKernel failed for VT prepass.\n");
            std::abort();
        }

        device->CopyBytesToHost(
            feedbackStatsHost, {feedbackStatsBuffer.data(), feedbackStatsBuffer.size()}, feedbackStatsBytes);
        const unsigned int sampledCount = feedbackStatsHost[0];
        const unsigned int overflowCount = feedbackStatsHost[1];
        const unsigned int copyCount = std::min(sampledCount, (unsigned int)feedbackCapacity);
        if (copyCount > 0)
        {
            device->CopyBytesToHost(feedbackKeysHost.data(),
                                    {feedbackKeysBuffer.data(), feedbackKeysBuffer.size()},
                                    size_t(copyCount) * sizeof(unsigned long long));
        }
        const std::vector<ybi::texture::VirtualTextureFeedbackEntry> histogram =
            WriteFeedbackFile("prepass.txt", -1, sampledCount, copyCount, overflowCount);
        std::printf("virtual-texture prepass feedback: sampled=%u stored=%u overflow=%u\n",
                    sampledCount,
                    copyCount,
                    overflowCount);
        ybi::texture::VirtualTextureUpdateStats prepassStats = {};
        if (!virtualTextureManager.ProcessFeedback(
                feedbackKeysHost.data(), copyCount, &prepassStats, &vtError))
        {
            fprintf(stderr, "virtual-texture prepass update failed: %s\n", vtError.c_str());
            std::abort();
        }
        LogProcessFeedbackProfile("vt_process_feedback", -1, prepassStats);
        std::printf("virtual-texture cache: unique=%u hits=%u misses=%u uploads=%u evictions=%u failed=%u\n",
                    prepassStats.uniqueCount,
                    prepassStats.hits,
                    prepassStats.misses,
                    prepassStats.uploads,
                    prepassStats.evictions,
                    prepassStats.failed);
        PrintVirtualTextureMemoryStats(
            setupMemoryStats, virtualTextureManager.GetMemoryStats(), renderDeviceBytes, renderHostBytes);
        virtualTextureManager.BindLaunchParams(&params);

        params.integrator = integrator == IntegratorType::AO
                                ? 1
                                : (integrator == IntegratorType::Path ? 3 : 0);
        params.feedbackSamplePercent = 10;
    }
    else
    {
        std::printf("memory: render resident\n");
        PrintByteStat("render.host", renderHostBytes);
        PrintByteStat("render.device", renderDeviceBytes);
    }

    for (int sppIndex = 0; sppIndex < sppPassCount; ++sppIndex)
    {
        const Clock::time_point passStart = Clock::now();
        params.spp = integrator == IntegratorType::AO ? 1 : 1;
        params.currentSpp = sppIndex;
        device->CopyBytesToDevice(paramsBuffer, &params, sizeof(LaunchParams));
        device->CopyBytesToDevice(feedbackStatsBuffer, feedbackStatsZero, feedbackStatsBytes);

        DispatchParams dispatchParams = {};
        dispatchParams.width = static_cast<uint32_t>(width);
        dispatchParams.height = static_cast<uint32_t>(height);
        dispatchParams.spp = static_cast<uint32_t>(params.spp);
        dispatchParams.launchParamsDevice = reinterpret_cast<uint64_t>(paramsBuffer.data());
        dispatchParams.launchParamsSize = sizeof(LaunchParams);
        dispatchParams.outputRGBA8 = imageBuffer;

        const RenderKernelId kernelId =
            integrator == IntegratorType::AO
                ? RenderKernelId::AO
                : (integrator == IntegratorType::Path ? RenderKernelId::PathTrace
                                                      : RenderKernelId::PrimaryDiffuse);
        if (!device->DispatchKernel(kernelId, dispatchParams))
        {
            fprintf(stderr, "DispatchKernel failed.\n");
            std::abort();
        }

        device->CopyBytesToHost(hostPassImage.data(), {imageBuffer.data(), imageBuffer.size()}, imageSize);
        const size_t pixelCount = static_cast<size_t>(width) * static_cast<size_t>(height);
        for (size_t i = 0; i < pixelCount; ++i)
        {
            accumRgb[i * 3 + 0] += hostPassImage[i * 4 + 0] / 255.0f;
            accumRgb[i * 3 + 1] += hostPassImage[i * 4 + 1] / 255.0f;
            accumRgb[i * 3 + 2] += hostPassImage[i * 4 + 2] / 255.0f;
        }

        device->CopyBytesToHost(
            feedbackStatsHost, {feedbackStatsBuffer.data(), feedbackStatsBuffer.size()}, feedbackStatsBytes);
        const unsigned int sampledCount = feedbackStatsHost[0];
        const unsigned int overflowCount = feedbackStatsHost[1];
        const unsigned int copyCount = std::min(sampledCount, (unsigned int)feedbackCapacity);
        if (copyCount > 0)
        {
            device->CopyBytesToHost(feedbackKeysHost.data(),
                                    {feedbackKeysBuffer.data(), feedbackKeysBuffer.size()},
                                    size_t(copyCount) * sizeof(unsigned long long));
        }

        char feedbackFileName[256];
        std::snprintf(feedbackFileName, sizeof(feedbackFileName), "spp_%04d.txt", sppIndex);
        (void)WriteFeedbackFile(
            feedbackFileName, sppIndex, sampledCount, copyCount, overflowCount);

        if (virtualTexture)
        {
            std::string vtError;
            ybi::texture::VirtualTextureUpdateStats passStats = {};
            if (!virtualTextureManager.ProcessFeedback(
                    feedbackKeysHost.data(), copyCount, &passStats, &vtError))
            {
                fprintf(stderr, "virtual-texture pass update failed: %s\n", vtError.c_str());
                std::abort();
            }
            LogProcessFeedbackProfile("vt_process_feedback", sppIndex, passStats);
            std::printf("virtual-texture cache spp=%d: unique=%u hits=%u misses=%u uploads=%u "
                        "evictions=%u failed=%u\n",
                        sppIndex,
                        passStats.uniqueCount,
                        passStats.hits,
                        passStats.misses,
                        passStats.uploads,
                        passStats.evictions,
                        passStats.failed);
            virtualTextureManager.BindLaunchParams(&params);
        }

        LogDuration("spp_pass_total",
                    sppIndex,
                    std::chrono::duration<double, std::milli>(Clock::now() - passStart).count());
    }

    std::vector<uint8_t> hostImage(imageSize, 0);
    const float invSpp = 1.0f / float(sppPassCount);
    for (size_t i = 0, pixelCount = static_cast<size_t>(width) * static_cast<size_t>(height);
         i < pixelCount;
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

    device->FreeBytes(paramsBuffer);
    virtualTextureManager.Shutdown();
    device->FreeBytes(feedbackStatsBuffer);
    device->FreeBytes(feedbackKeysBuffer);
    device->FreeBytes(imageBuffer);
    printf("render: end\n");
    fflush(stdout);
}
} // namespace

struct HarnessState
{
    ScenePool scenePool = {};
    ScenePool flattenedScenePool = {};
    Scene *rootScene = nullptr;

    DeviceMemoryView<uint8_t> materialsBuffer = {};
    DeviceMemoryView<uint8_t> lightsBuffer = {};
    DeviceMemoryView<uint8_t> lightShadowExcludeRefsBuffer = {};
    DeviceMemoryView<uint8_t> materialTextureRefsBuffer = {};
    std::vector<PackedMaterial> packedMaterials;
    std::vector<PackedLight> packedLights;
    std::vector<int> lightShadowExcludeRefs;
    int materialTextureRefCount = 0;
    int materialTextureRefStride = kUdimSlotCount;
    int materialTextureRefSemanticCount = kMaterialSemanticCount;
    std::vector<LaunchParams::MaterialTextureRef> materialTextureRefsHost;
    std::vector<DecodedMaterialTexture> decodedTextures;
    UploadedMaterialTextures uploadedMaterialTextures = {};
    DecodedMaterialTexture domeTexture = {};
    UploadedStandaloneTexture uploadedDomeTexture = {};
    std::vector<ybi::texture::VirtualTextureRegisterInput> virtualTextureRegistrations;

    std::vector<SceneMeshUploadRef> meshUploadRefs;
    UploadedMeshRefs uploadedRefs = {};
    DeviceMemoryView<uint8_t> instanceGeomRefsBuffer = {};
    HarnessMemoryStats memoryStats = {};

    std::optional<RenderCameraOverride> usdCamera = std::nullopt;
};

static bool UploadScenePhase(Device *device, const CliOptions &options, HarnessState *state)
{
    YBI_ASSERT(device);
    YBI_ASSERT(state);

    using Clock = std::chrono::steady_clock;
    const auto phaseStart = Clock::now();
    auto LogPhase = [](const char *label, Clock::time_point start) {
        const auto now = Clock::now();
        const double ms =
            std::chrono::duration<double, std::milli>(now - start).count();
        std::printf("profile: %s %.3f ms\n", label, ms);
        return now;
    };

    printf("template_harness: loading usd scene %s\n", options.inputPath.c_str());
    fflush(stdout);

    USDLoadOptions loadOptions = {};
    loadOptions.purposes = options.purposes;
    loadOptions.camera = options.usdCamera;
    LoadUSDScene(&state->scenePool, options.inputPath, loadOptions);
    auto tLoad = LogPhase("usd_load", phaseStart);
    if (state->scenePool.scenes.empty() ||
        state->scenePool.rootSceneIndex >= state->scenePool.scenes.size())
    {
        fprintf(
            stderr, "Failed to load USD scene or invalid root: %s\n", options.inputPath.c_str());
        return false;
    }

    std::string flattenError;
    if (!FlattenScenePoolToRootChildren(
            &state->scenePool, &state->flattenedScenePool, &flattenError))
    {
        fprintf(stderr, "Failed to flatten USD ScenePool: %s\n", flattenError.c_str());
        return false;
    }
    auto tFlatten = LogPhase("flatten_scene_pool", tLoad);
    if (state->flattenedScenePool.scenes.empty() ||
        state->flattenedScenePool.rootSceneIndex >= state->flattenedScenePool.scenes.size())
    {
        fprintf(
            stderr, "Flattened ScenePool invalid for USD scene: %s\n", options.inputPath.c_str());
        return false;
    }
    state->rootScene =
        state->flattenedScenePool.scenes[state->flattenedScenePool.rootSceneIndex].get();
    YBI_ASSERT(state->rootScene);
    if (!TessellateRootSubdivisionMeshes(state->rootScene, state->flattenedScenePool.camera))
    {
        return false;
    }
    auto tTess = LogPhase("tessellate_subdivision", tFlatten);

    state->memoryStats = {};
    state->memoryStats.hostSourceSceneBytes = ComputeScenePoolHostBytes(state->scenePool);
    state->memoryStats.hostFlattenedSceneBytes = ComputeScenePoolHostBytes(state->flattenedScenePool);
    state->memoryStats.hostTessellatedMeshBytes = ComputeMeshVectorHostBytes(state->rootScene->meshes);

    DestroyUploadedTexture(device, &state->uploadedDomeTexture);
    state->decodedTextures.clear();
    state->domeTexture = {};
    state->packedMaterials.clear();
    state->packedLights.clear();
    state->lightShadowExcludeRefs.clear();
    state->materialTextureRefsHost.clear();
    state->virtualTextureRegistrations.clear();

    state->packedMaterials.reserve(state->scenePool.materials.size());
    for (const MaterialInfo &material : state->scenePool.materials)
    {
        state->packedMaterials.push_back(material.packed);
    }
    if (!state->packedMaterials.empty())
    {
        const size_t bytes = state->packedMaterials.size() * sizeof(PackedMaterial);
        state->materialsBuffer = device->AllocBytes(bytes);
        state->memoryStats.deviceMaterialBytes = state->materialsBuffer.numBytes();
        device->CopyBytesToDevice(state->materialsBuffer, state->packedMaterials.data(), bytes);
    }

    std::vector<LaunchParams::MaterialTextureRef> launchRefs(
        state->scenePool.materials.size() * static_cast<size_t>(kMaterialSemanticCount) *
        static_cast<size_t>(kUdimSlotCount));
    Clock::time_point tTextureReady = tTess;

    if (options.virtualTexture)
    {
        if (options.useNtc)
        {
            std::fprintf(stderr,
                         "virtual-texture: ignoring --ntc; strict VT path uses tile bins only\n");
        }
        if (options.virtualTextureTilesDir.empty())
        {
            std::fprintf(stderr, "virtual-texture: --vt-tiles-dir is required\n");
            return false;
        }

        std::vector<std::pair<unsigned int, std::string>> texturePaths;
        std::vector<ybi::texture::VirtualTextureMaterialSource> vtSources;
        texturePaths.reserve(state->scenePool.materials.size() * static_cast<size_t>(kMaterialSemanticCount));
        vtSources.reserve(state->scenePool.materials.size() * static_cast<size_t>(kMaterialSemanticCount));
        for (size_t materialIndex = 0; materialIndex < state->scenePool.materials.size(); ++materialIndex)
        {
            for (const MaterialTextureInput &input : state->scenePool.materials[materialIndex].textureInputs)
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
                const unsigned int textureId = static_cast<unsigned int>(
                    materialIndex * static_cast<size_t>(kMaterialSemanticCount) +
                    static_cast<size_t>(semantic));
                texturePaths.emplace_back(textureId, input.texturePath);
                ybi::texture::VirtualTextureMaterialSource source = {};
                source.textureId = textureId;
                source.materialIndex = static_cast<uint32_t>(materialIndex);
                source.semanticIndex = static_cast<int>(semantic);
                source.texturePath = input.texturePath;
                source.wrapS = input.wrapS;
                source.wrapT = input.wrapT;
                vtSources.push_back(std::move(source));
            }
        }

        const std::unordered_map<unsigned int, std::string> tileFiles =
            ybi::texture::BuildVirtualTextureTileFileMap(texturePaths, options.virtualTextureTilesDir);
        ybi::texture::VirtualTextureMaterialBuildResult vtBuild = {};
        std::string vtError;
        if (!ybi::texture::BuildVirtualTextureMaterialMetadata(state->scenePool.materials.size(),
                                                               state->materialTextureRefSemanticCount,
                                                               state->materialTextureRefStride,
                                                               vtSources,
                                                               tileFiles,
                                                               &vtBuild,
                                                               &vtError))
        {
            std::fprintf(stderr, "%s\n", vtError.c_str());
            return false;
        }
        launchRefs = std::move(vtBuild.materialTextureRefs);
        state->virtualTextureRegistrations = std::move(vtBuild.registrations);
        state->memoryStats.virtualTextureVirtualBytes = vtBuild.totalVirtualTextureBytes;
        std::printf("virtual-texture: mapped %u texture slots and %u active UDIMs from %s\n",
                    vtBuild.mappedMaterialCount,
                    vtBuild.activeUdimCount,
                    options.virtualTextureTilesDir.c_str());
        if (vtBuild.skippedMissingTileCount > 0u)
        {
            std::printf("virtual-texture: skipped %u missing tile bins\n",
                        vtBuild.skippedMissingTileCount);
        }
        std::printf(
            "virtual-texture: strict VT path enabled; skipping image decode and classic texture upload\n");
        tTextureReady = LogPhase("build_virtual_texture_metadata", tTess);
    }
    else
    {
        if (!DecodeImageTextures(state->scenePool.materials, &state->decodedTextures))
        {
            fprintf(stderr, "Image runtime decode failed.\n");
            return false;
        }
#if defined(YBI_OPTIX_HARNESS_WITH_NTC)
        if (options.useNtc)
        {
            std::vector<testbvh::DecodedMaterialTexture> ntcTextures;
            std::string ntcError;
            if (testbvh::DecodeNtcDiffuseTextures(state->scenePool.materials, &ntcTextures, &ntcError))
            {
                int overrideCount = 0;
                const size_t numMaterials =
                    std::min(state->scenePool.materials.size(), ntcTextures.size());
                for (size_t i = 0; i < numMaterials; ++i)
                {
                    if (!ntcTextures[i].valid)
                    {
                        continue;
                    }
                    const int slot =
                        MaterialTextureSlotIndex(i, MaterialTextureSemantic::Diffuse, kUdimMin);
                    YBI_ASSERT(slot >= 0 && static_cast<size_t>(slot) < state->decodedTextures.size());
                    DecodedMaterialTexture &dst = state->decodedTextures[static_cast<size_t>(slot)];
                    dst.valid = true;
                    dst.udim = kUdimMin;
                    dst.width = ntcTextures[i].width;
                    dst.height = ntcTextures[i].height;
                    dst.format = DeviceTextureFormat::RGBA8_UNORM;
                    dst.pixels = std::move(ntcTextures[i].rgba8);
                    dst.sourcePath = ntcTextures[i].ntcPath;
                    const MaterialTextureInput *diffuse = FindTextureInputBySemantic(
                        state->scenePool.materials[i], MaterialTextureSemantic::Diffuse);
                    dst.wrapS = diffuse ? diffuse->wrapS : TEXTURE_WRAP_MODE_REPEAT;
                    dst.wrapT = diffuse ? diffuse->wrapT : TEXTURE_WRAP_MODE_REPEAT;
                    overrideCount++;
                }
                printf("NTC runtime: applied %d material overrides\n", overrideCount);
            }
            else
            {
                fprintf(stderr,
                        "NTC runtime decode failed (continuing with image textures): %s\n",
                        ntcError.c_str());
            }
        }
#else
        if (options.useNtc)
        {
            fprintf(stderr,
                    "NTC runtime requested via --ntc, but harness built without WITH_NTC. "
                    "Using image textures only.\n");
        }
#endif
        state->memoryStats.hostDecodedTextureBytes = 0u;
        for (const DecodedMaterialTexture &texture : state->decodedTextures)
        {
            state->memoryStats.hostDecodedTextureBytes += static_cast<uint64_t>(texture.pixels.size());
        }

        tTextureReady = LogPhase("decode_image_textures", tTess);
        std::string textureUploadError;
        if (!UploadDecodedTextures(
                device, state->decodedTextures, &state->uploadedMaterialTextures, &textureUploadError))
        {
            fprintf(stderr, "Texture upload failed: %s\n", textureUploadError.c_str());
            return false;
        }
        state->memoryStats.deviceUploadedTextureBytes = state->uploadedMaterialTextures.deviceBytes;
        LogPhase("upload_textures", tTextureReady);
        for (size_t i = 0; i < state->uploadedMaterialTextures.refs.size(); ++i)
        {
            const LaunchParams::MaterialTextureRef &src = state->uploadedMaterialTextures.refs[i];
            launchRefs[i] = src;
        }
    }

    if (!launchRefs.empty())
    {
        const size_t refsBytes = launchRefs.size() * sizeof(LaunchParams::MaterialTextureRef);
        state->materialTextureRefsHost = launchRefs;
        state->memoryStats.hostMaterialTextureRefBytes = refsBytes;
        state->materialTextureRefsBuffer = device->AllocBytes(refsBytes);
        state->memoryStats.deviceMaterialTextureRefBytes = state->materialTextureRefsBuffer.numBytes();
        device->CopyBytesToDevice(state->materialTextureRefsBuffer, launchRefs.data(), refsBytes);
        state->materialTextureRefCount = static_cast<int>(state->scenePool.materials.size());
    }

    std::string domeTexturePath;
    if (FindFirstDomeTexturePath(state->scenePool.lights, &domeTexturePath))
    {
        state->domeTexture.valid = true;
        state->domeTexture.udim = kUdimMin;
        state->domeTexture.sourcePath = domeTexturePath;
        state->domeTexture.wrapS = TEXTURE_WRAP_MODE_REPEAT;
        state->domeTexture.wrapT = TEXTURE_WRAP_MODE_CLAMP;
        std::string domeReason;
        const bool domeIsExr = ybi::texture::LowerExt(domeTexturePath) == ".exr";
        state->domeTexture.format =
            domeIsExr ? DeviceTextureFormat::RGBA16_FLOAT : DeviceTextureFormat::RGBA8_UNORM;
        const bool domeLoaded =
            domeIsExr
                ? LoadExrRgba16Float(domeTexturePath,
                                     &state->domeTexture.width,
                                     &state->domeTexture.height,
                                     &state->domeTexture.pixels,
                                     &domeReason,
                                     false)
                : LoadImageRgba8(domeTexturePath,
                                 &state->domeTexture.width,
                                 &state->domeTexture.height,
                                 &state->domeTexture.pixels,
                                 &domeReason,
                                 false);
        if (!domeLoaded)
        {
            std::fprintf(stderr, "Dome texture load failed (%s): %s\n", domeTexturePath.c_str(), domeReason.c_str());
            return false;
        }

        state->memoryStats.hostDecodedTextureBytes +=
            static_cast<uint64_t>(state->domeTexture.pixels.size());
        std::string domeUploadError;
        if (!UploadDecodedTexture(
                device, state->domeTexture, &state->uploadedDomeTexture, &domeUploadError))
        {
            std::fprintf(stderr, "Dome texture upload failed: %s\n", domeUploadError.c_str());
            return false;
        }
        state->memoryStats.deviceUploadedTextureBytes += state->uploadedDomeTexture.deviceBytes;
    }
    auto tTexRefs = LogPhase("upload_texture_refs", tTextureReady);

    CollectScenePoolMeshUploadRefs(&state->flattenedScenePool, &state->meshUploadRefs);
    auto tCollect = LogPhase("collect_mesh_upload_refs", tTexRefs);
    ResolveLightShadowExcludes(
        state->meshUploadRefs, state->scenePool.lights, &state->packedLights, &state->lightShadowExcludeRefs);
    for (size_t i = 0; i < state->scenePool.lights.size() && i < state->packedLights.size(); ++i)
    {
        if (state->packedLights[i].shadowExcludeCount == 0u)
        {
            continue;
        }
        std::printf("light shadow excludes: %s -> %u mesh refs\n",
                    state->scenePool.lights[i].lightPath.c_str(),
                    state->packedLights[i].shadowExcludeCount);
    }
    if (!state->packedLights.empty())
    {
        const size_t bytes = state->packedLights.size() * sizeof(PackedLight);
        state->lightsBuffer = device->AllocBytes(bytes);
        state->memoryStats.deviceLightBytes = state->lightsBuffer.numBytes();
        device->CopyBytesToDevice(state->lightsBuffer, state->packedLights.data(), bytes);
    }
    if (!state->lightShadowExcludeRefs.empty())
    {
        const size_t bytes = state->lightShadowExcludeRefs.size() * sizeof(int);
        state->lightShadowExcludeRefsBuffer = device->AllocBytes(bytes);
        state->memoryStats.deviceLightBytes += state->lightShadowExcludeRefsBuffer.numBytes();
        device->CopyBytesToDevice(
            state->lightShadowExcludeRefsBuffer, state->lightShadowExcludeRefs.data(), bytes);
    }
    auto bvhStart = Clock::now();
    for (const std::unique_ptr<Scene> &scenePtr : state->flattenedScenePool.scenes)
    {
        Scene *scene = scenePtr.get();
        YBI_ASSERT(scene);
        device->BuildBVH(scene);
    }
    LogPhase("build_bvh_total", bvhStart);
    if (state->rootScene->bvhHandle == 0)
    {
        fprintf(stderr, "USD root scene BVH is invalid: %s\n", options.inputPath.c_str());
        return false;
    }
    auto tBvh = LogPhase("build_bvh_phase", tCollect);

    state->uploadedRefs = UploadScenePoolMeshRefs(device, state->meshUploadRefs);
    state->memoryStats.deviceMeshUploadBytes = state->uploadedRefs.deviceBytes;
    auto tUploadRefs = LogPhase("upload_scene_mesh_refs", tBvh);
    if (!state->uploadedRefs.refs.empty())
    {
        const size_t refsBytes =
            state->uploadedRefs.refs.size() * sizeof(LaunchParams::InstanceGeomRef);
        state->instanceGeomRefsBuffer = device->AllocBytes(refsBytes);
        state->memoryStats.deviceInstanceRefBytes = state->instanceGeomRefsBuffer.numBytes();
        device->CopyBytesToDevice(
            state->instanceGeomRefsBuffer, state->uploadedRefs.refs.data(), refsBytes);
    }
    auto tInstanceRefs = LogPhase("upload_instance_geom_refs", tUploadRefs);
    state->memoryStats.deviceBvhBytes = device->GetBVHAllocatedBytes();

    if (!options.cameraPosition.has_value() && !options.lookAt.has_value())
    {
        state->usdCamera = BuildUsdRenderCamera(state->flattenedScenePool.camera);
        if (!state->usdCamera.has_value())
        {
            fprintf(stderr, "USD camera missing/invalid: %s\n", options.inputPath.c_str());
            return false;
        }
        printf("template_harness: using usd camera %s viewport=%dx%d\n",
               state->flattenedScenePool.camera.path.c_str(),
               state->usdCamera->width,
               state->usdCamera->height);
        fflush(stdout);
    }

    PrintSetupMemoryStats(state->memoryStats);
    LogPhase("upload_scene_total", phaseStart);
    return true;
}

static bool RenderPhase(Device *device, const CliOptions &options, const HarnessState &state)
{
    YBI_ASSERT(device);
    if (device->GetKind() == DeviceKind::GPU)
    {
        const std::string ptx = ReadTextFile(YBI_OPTIX_PRIMARY_PTX_PATH);
        printf("template_harness: ptx loaded\n");
        fflush(stdout);

        std::string kernelInitError;
        if (!device->InitializeKernels(ptx, &kernelInitError))
        {
            fprintf(stderr, "Failed to initialize kernels: %s\n", kernelInitError.c_str());
            return false;
        }
        printf("template_harness: pipeline created\n");
        fflush(stdout);
    }
    else
    {
        std::string kernelInitError;
        if (!device->InitializeKernels("", &kernelInitError))
        {
            fprintf(stderr, "Failed to initialize CPU kernels: %s\n", kernelInitError.c_str());
            return false;
        }
    }

    const ybi::Vec3 dummyBoundsMin = ybi::Vec3(-1.0f, -1.0f, -1.0f);
    const ybi::Vec3 dummyBoundsMax = ybi::Vec3(1.0f, 1.0f, 1.0f);
    RenderTraversable(device,
                      state.rootScene->bvhHandle,
                      dummyBoundsMin,
                      dummyBoundsMax,
                      options.outputPath.c_str(),
                      options.integrator,
                      options.spp,
                      options.maxDepth,
                      reinterpret_cast<DevicePtr>(state.instanceGeomRefsBuffer.data()),
                      (int)state.uploadedRefs.refs.size(),
                      reinterpret_cast<DevicePtr>(state.materialsBuffer.data()),
                      static_cast<int>(state.packedMaterials.size()),
                      reinterpret_cast<DevicePtr>(state.lightsBuffer.data()),
                      static_cast<int>(state.packedLights.size()),
                      reinterpret_cast<DevicePtr>(state.lightShadowExcludeRefsBuffer.data()),
                      static_cast<int>(state.lightShadowExcludeRefs.size()),
                      reinterpret_cast<DevicePtr>(state.materialTextureRefsBuffer.data()),
                      state.materialTextureRefCount,
                      state.materialTextureRefStride,
                      state.materialTextureRefSemanticCount,
                      state.uploadedDomeTexture.ref,
                      options.textureView,
                      options.virtualTexture,
                      state.virtualTextureRegistrations,
                      state.memoryStats,
                      options.virtualTextureCacheBytes,
                      options.virtualTextureTailMaxDim,
                      options.virtualTextureMaxPageUploads,
                      options.singlePixelMode,
                      options.singlePixelX,
                      options.singlePixelY,
                      options.writeFeedbackFiles,
                      state.scenePool.camera.upAxis,
                      state.usdCamera,
                      options.cameraPosition,
                      options.lookAt);
    printf("Wrote %s\n", options.outputPath.c_str());
    return true;
}

static void DeinitPhase(Device *device, HostMemoryArena *hostArena, HarnessState *state)
{
    YBI_ASSERT(device);
    YBI_ASSERT(hostArena);
    YBI_ASSERT(state);

    if (state->instanceGeomRefsBuffer.data())
    {
        device->FreeBytes(state->instanceGeomRefsBuffer);
    }
    if (state->materialsBuffer.data())
    {
        device->FreeBytes(state->materialsBuffer);
    }
    if (state->lightsBuffer.data())
    {
        device->FreeBytes(state->lightsBuffer);
    }
    if (state->lightShadowExcludeRefsBuffer.data())
    {
        device->FreeBytes(state->lightShadowExcludeRefsBuffer);
    }
    if (state->materialTextureRefsBuffer.data())
    {
        device->FreeBytes(state->materialTextureRefsBuffer);
    }
    DestroyUploadedTexture(device, &state->uploadedDomeTexture);
    DestroyUploadedTextures(device, &state->uploadedMaterialTextures);
    for (DeviceMemoryView<uint8_t> &buffer : state->uploadedRefs.ownedBuffers)
    {
        device->FreeBytes(buffer);
    }
    hostArena->Clear();
    device->ClearTransientMemory();
    device->DestroyKernels();
}

int main(int argc, char **argv)
{
    printf("template_harness: start\n");
    fflush(stdout);
    const CliOptions options = ParseCli(argc, argv);
    printf("template_harness: parsed cli\n");
    fflush(stdout);

    std::unique_ptr<Device> deviceStorage;
    Device *device = Device::CreateDevice(options.deviceKind, deviceStorage);
    if (!device)
    {
        fprintf(stderr, "Failed to create device.\n");
        return 1;
    }
    if (device->GetKind() == DeviceKind::GPU)
    {
        printf("template_harness: cuda device ready\n");
    }
    else
    {
        printf("template_harness: cpu device ready\n");
    }
    fflush(stdout);
    HostMemoryArena hostArena;
    HarnessState state = {};
    const bool uploadOk = UploadScenePhase(device, options, &state);
    bool renderOk = false;
    if (uploadOk)
    {
        renderOk = RenderPhase(device, options, state);
    }
    DeinitPhase(device, &hostArena, &state);
    if (!uploadOk || !renderOk)
    {
        return 1;
    }
    std::fflush(stdout);
    std::fflush(stderr);
    std::_Exit(0);
}
