#include "device/cpu_device.h"
#include "io/usd/load.h"
#include "scene/scene.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "third_party/stb_image_write.h"
#include "util/float3.h"
#include "util/float3x4.h"
#include "util/float4.h"
#include "util/float4x4.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <optional>
#include <string>
#include <vector>

#if __has_include(<embree4/rtcore.h>)
#include <embree4/rtcore.h>
#elif __has_include(<embree3/rtcore.h>)
#include <embree3/rtcore.h>
#else
#error "Embree headers not found: expected embree3/rtcore.h or embree4/rtcore.h"
#endif

using namespace ybi;

namespace
{
enum class IntegratorType
{
    Primary,
    AO
};

struct CliOptions
{
    std::string inputPath;
    std::string outputPath = "embree_usd_scene.png";
    IntegratorType integrator = IntegratorType::AO;
    int spp = 4;
    float aoMaxDistance = 0.0f;
};

struct RenderCamera
{
    int width = 1280;
    int height = 720;
    float3 origin = make_float3(0.0f);
    float3 U = make_float3(1.0f, 0.0f, 0.0f);
    float3 V = make_float3(0.0f, 1.0f, 0.0f);
    float3 W = make_float3(0.0f, 0.0f, 1.0f);
};

struct Bounds
{
    float3 min = make_float3(std::numeric_limits<float>::max());
    float3 max = make_float3(-std::numeric_limits<float>::max());
    bool valid = false;
};

static void PrintUsage(const char *exeName)
{
    printf("Usage: %s --file <scene.usd[a|c]> [--out path.png] [--integrator primary|ao] "
           "[--spp N] [--ao-distance D]\n",
           exeName);
}

static std::string Lowercase(std::string value)
{
    std::transform(value.begin(), value.end(), value.begin(), [](unsigned char c) {
        return static_cast<char>(std::tolower(c));
    });
    return value;
}

static CliOptions ParseCli(int argc, char **argv)
{
    CliOptions options = {};
    for (int i = 1; i < argc; i++)
    {
        const std::string arg = argv[i];
        if (arg == "--file")
        {
            if (i + 1 >= argc)
            {
                PrintUsage(argv[0]);
                std::exit(1);
            }
            options.inputPath = argv[++i];
            continue;
        }
        if (arg == "--out")
        {
            if (i + 1 >= argc)
            {
                PrintUsage(argv[0]);
                std::exit(1);
            }
            options.outputPath = argv[++i];
            continue;
        }
        if (arg == "--integrator")
        {
            if (i + 1 >= argc)
            {
                PrintUsage(argv[0]);
                std::exit(1);
            }
            const std::string mode = argv[++i];
            if (mode == "primary")
            {
                options.integrator = IntegratorType::Primary;
            }
            else if (mode == "ao")
            {
                options.integrator = IntegratorType::AO;
            }
            else
            {
                PrintUsage(argv[0]);
                std::exit(1);
            }
            continue;
        }
        if (arg == "--spp")
        {
            if (i + 1 >= argc)
            {
                PrintUsage(argv[0]);
                std::exit(1);
            }
            options.spp = std::max(1, std::stoi(argv[++i]));
            continue;
        }
        if (arg == "--ao-distance")
        {
            if (i + 1 >= argc)
            {
                PrintUsage(argv[0]);
                std::exit(1);
            }
            options.aoMaxDistance = std::max(0.0f, std::stof(argv[++i]));
            continue;
        }
        if (arg == "--help" || arg == "-h")
        {
            PrintUsage(argv[0]);
            std::exit(0);
        }
        PrintUsage(argv[0]);
        std::exit(1);
    }

    if (options.inputPath.empty())
    {
        PrintUsage(argv[0]);
        std::exit(1);
    }

    const std::string lowered = Lowercase(options.inputPath);
    const bool isUsd = lowered.size() >= 4 &&
                       (lowered.rfind(".usd") == lowered.size() - 4 ||
                        lowered.rfind(".usda") == lowered.size() - 5 ||
                        lowered.rfind(".usdc") == lowered.size() - 5);
    if (!isUsd)
    {
        fprintf(stderr, "Embree harness only supports USD files: %s\n", options.inputPath.c_str());
        std::exit(1);
    }
    if (options.integrator == IntegratorType::Primary)
    {
        options.spp = 1;
    }
    return options;
}

static float3 Cross(const float3 &a, const float3 &b)
{
    return make_float3(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
}

static float3 Normalize(const float3 &v)
{
    const float lenSq = dot(v, v);
    if (lenSq <= 1e-16f)
    {
        return make_float3(0.0f, 0.0f, 1.0f);
    }
    return v * (1.0f / std::sqrt(lenSq));
}

static float4x4 ToFloat4x4(const float3x4 &m)
{
    return float4x4(m.m[0][0],
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

static bool InvertAffine(const float4x4 &m, float4x4 &out)
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

    out = float4x4(
        r00, r01, r02, itx, r10, r11, r12, ity, r20, r21, r22, itz, 0.0f, 0.0f, 0.0f, 1.0f);
    return true;
}

static std::optional<RenderCamera> BuildUsdRenderCamera(const Camera &camera)
{
    if (camera.viewportWidth <= 0 || camera.viewportHeight <= 0)
    {
        return std::nullopt;
    }

    float4x4 worldFromCamera = {};
    if (!InvertAffine(camera.cameraFromWorld, worldFromCamera))
    {
        return std::nullopt;
    }

    float3 right = Normalize(
        make_float3(worldFromCamera.m[0][0], worldFromCamera.m[1][0], worldFromCamera.m[2][0]));
    float3 up = Normalize(
        make_float3(worldFromCamera.m[0][1], worldFromCamera.m[1][1], worldFromCamera.m[2][1]));
    const float3 forward =
        Normalize(make_float3(-worldFromCamera.m[0][2],
                              -worldFromCamera.m[1][2],
                              -worldFromCamera.m[2][2]));

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

    RenderCamera out = {};
    out.width = camera.viewportWidth;
    out.height = camera.viewportHeight;
    out.origin =
        make_float3(worldFromCamera.m[0][3], worldFromCamera.m[1][3], worldFromCamera.m[2][3]);
    out.U = right * (aspect * tanHalfFov);
    out.V = up * tanHalfFov;
    out.W = forward;
    return out;
}

static void ExpandBounds(Bounds &bounds, const float3 &p)
{
    bounds.min.x = std::min(bounds.min.x, p.x);
    bounds.min.y = std::min(bounds.min.y, p.y);
    bounds.min.z = std::min(bounds.min.z, p.z);
    bounds.max.x = std::max(bounds.max.x, p.x);
    bounds.max.y = std::max(bounds.max.y, p.y);
    bounds.max.z = std::max(bounds.max.z, p.z);
    bounds.valid = true;
}

static void ComputeWorldBounds(ScenePool &scenePool, Bounds &outBounds)
{
    outBounds = {};
    if (scenePool.scenes.empty() || scenePool.rootSceneIndex >= scenePool.scenes.size())
    {
        return;
    }

    struct StackEntry
    {
        Scene *scene = nullptr;
        float4x4 worldFromScene = float4x4::Identity();
    };

    std::vector<StackEntry> stack;
    stack.push_back({scenePool.scenes[scenePool.rootSceneIndex].get(), float4x4::Identity()});
    while (!stack.empty())
    {
        const StackEntry entry = stack.back();
        stack.pop_back();
        if (!entry.scene)
        {
            continue;
        }

        for (const Mesh &mesh : entry.scene->meshes)
        {
            const float4x4 worldFromMesh = mul(entry.worldFromScene, ToFloat4x4(mesh.parentFromLocal));
            for (const float3 &p : mesh.positions)
            {
                const float4 wp = mul(worldFromMesh, make_float4(p.x, p.y, p.z, 1.0f));
                ExpandBounds(outBounds, make_float3(wp.x, wp.y, wp.z));
            }
        }

        for (const Instance &instance : entry.scene->instances)
        {
            if (instance.childSceneIndex >= entry.scene->childScenes.size())
            {
                continue;
            }
            Scene *child = entry.scene->childScenes[instance.childSceneIndex];
            stack.push_back({child, mul(entry.worldFromScene, ToFloat4x4(instance.parentFromLocal))});
        }
    }
}

static uint32_t Hash32(uint32_t x)
{
    x ^= x >> 17;
    x *= 0xed5ad4bbU;
    x ^= x >> 11;
    x *= 0xac4c1b51U;
    x ^= x >> 15;
    x *= 0x31848babU;
    x ^= x >> 14;
    return x;
}

static float3 MakeFlatAlbedo(uint32_t id)
{
    const uint32_t h = Hash32(id + 1u);
    const float r = 0.2f + 0.8f * ((h & 0xFFu) / 255.0f);
    const float g = 0.2f + 0.8f * (((h >> 8) & 0xFFu) / 255.0f);
    const float b = 0.2f + 0.8f * (((h >> 16) & 0xFFu) / 255.0f);
    return make_float3(r, g, b);
}

static float NextFloat01(uint32_t &state)
{
    state = 1664525u * state + 1013904223u;
    return static_cast<float>(state & 0x00FFFFFFu) / static_cast<float>(0x01000000u);
}

static float3 SampleCosineHemisphere(const float3 &n, uint32_t &rngState)
{
    const float r1 = NextFloat01(rngState);
    const float r2 = NextFloat01(rngState);
    const float phi = 2.0f * 3.14159265358979323846f * r1;
    const float r = std::sqrt(r2);
    const float x = r * std::cos(phi);
    const float y = r * std::sin(phi);
    const float z = std::sqrt(std::max(0.0f, 1.0f - r2));

    const float3 up = std::abs(n.z) < 0.999f ? make_float3(0.0f, 0.0f, 1.0f)
                                             : make_float3(0.0f, 1.0f, 0.0f);
    const float3 tangent = Normalize(Cross(up, n));
    const float3 bitangent = Cross(n, tangent);
    return Normalize(tangent * x + bitangent * y + n * z);
}

static bool SavePNG(const char *filePath, const std::vector<uint8_t> &rgba, int width, int height)
{
    const int strideInBytes = width * 4;
    return stbi_write_png(filePath, width, height, 4, rgba.data(), strideInBytes) != 0;
}

static float Clamp01(float x)
{
    if (x < 0.0f)
    {
        return 0.0f;
    }
    if (x > 1.0f)
    {
        return 1.0f;
    }
    return x;
}

static std::vector<uint8_t> RenderImage(RTCScene rootScene,
                                        const RenderCamera &camera,
                                        IntegratorType integrator,
                                        int spp,
                                        float aoBias,
                                        float aoMaxDistance,
                                        const std::vector<float3> &albedoByRefIndex)
{
    const int width = camera.width;
    const int height = camera.height;
    std::vector<uint8_t> image(static_cast<size_t>(width) * static_cast<size_t>(height) * 4u, 0u);
    for (int y = 0; y < height; y++)
    {
        for (int x = 0; x < width; x++)
        {
            float3 pixelColor = make_float3(0.0f);
            const int sampleCount = integrator == IntegratorType::AO ? std::max(1, spp) : 1;
            for (int s = 0; s < sampleCount; s++)
            {
                uint32_t rngState = Hash32(static_cast<uint32_t>((y * width + x) * 9781 + s * 6271 + 1));
                const float jitterX = integrator == IntegratorType::AO ? NextFloat01(rngState) : 0.5f;
                const float jitterY = integrator == IntegratorType::AO ? NextFloat01(rngState) : 0.5f;
                const float u = ((static_cast<float>(x) + jitterX) / static_cast<float>(width)) * 2.0f -
                                1.0f;
                const float v = 1.0f -
                                ((static_cast<float>(y) + jitterY) / static_cast<float>(height)) * 2.0f;
                const float3 rayDir = Normalize(camera.W + camera.U * u + camera.V * v);

                RTCRayHit rayHit = {};
                rayHit.ray.org_x = camera.origin.x;
                rayHit.ray.org_y = camera.origin.y;
                rayHit.ray.org_z = camera.origin.z;
                rayHit.ray.dir_x = rayDir.x;
                rayHit.ray.dir_y = rayDir.y;
                rayHit.ray.dir_z = rayDir.z;
                rayHit.ray.tnear = 0.0f;
                rayHit.ray.tfar = std::numeric_limits<float>::infinity();
                rayHit.ray.mask = ~0u;
                rayHit.ray.flags = 0u;
                rayHit.hit.geomID = RTC_INVALID_GEOMETRY_ID;
#if YBI_EMBREE_VERSION_MAJOR >= 4
                RTCRayQueryContext rayQueryContext;
                rtcInitRayQueryContext(&rayQueryContext);
                RTCIntersectArguments intersectArgs;
                rtcInitIntersectArguments(&intersectArgs);
                intersectArgs.context = &rayQueryContext;
                rtcIntersect1(rootScene, &rayHit, &intersectArgs);
#else
                RTCIntersectContext context;
                rtcInitIntersectContext(&context);
                rtcIntersect1(rootScene, &context, &rayHit);
#endif
                if (rayHit.hit.geomID == RTC_INVALID_GEOMETRY_ID)
                {
                    pixelColor += make_float3(0.85f, 0.85f, 0.9f);
                    continue;
                }

                const uint32_t refIndex = rayHit.hit.geomID;
                float3 albedo = make_float3(0.7f);
                if (refIndex < albedoByRefIndex.size())
                {
                    albedo = albedoByRefIndex[refIndex];
                }

                if (integrator == IntegratorType::Primary)
                {
                    pixelColor += albedo;
                    continue;
                }

                const float3 n = Normalize(make_float3(rayHit.hit.Ng_x, rayHit.hit.Ng_y, rayHit.hit.Ng_z));
                const float3 hitPoint = camera.origin + rayDir * rayHit.ray.tfar;
                const float3 aoDir = SampleCosineHemisphere(n, rngState);

                RTCRay shadowRay = {};
                shadowRay.org_x = hitPoint.x + n.x * aoBias;
                shadowRay.org_y = hitPoint.y + n.y * aoBias;
                shadowRay.org_z = hitPoint.z + n.z * aoBias;
                shadowRay.dir_x = aoDir.x;
                shadowRay.dir_y = aoDir.y;
                shadowRay.dir_z = aoDir.z;
                shadowRay.tnear = 0.0f;
                shadowRay.tfar = aoMaxDistance;
                shadowRay.mask = ~0u;
                shadowRay.flags = 0u;

#if YBI_EMBREE_VERSION_MAJOR >= 4
                RTCRayQueryContext occlusionQueryContext;
                rtcInitRayQueryContext(&occlusionQueryContext);
                RTCOccludedArguments occludedArgs;
                rtcInitOccludedArguments(&occludedArgs);
                occludedArgs.context = &occlusionQueryContext;
                rtcOccluded1(rootScene, &shadowRay, &occludedArgs);
#else
                RTCIntersectContext shadowContext;
                rtcInitIntersectContext(&shadowContext);
                rtcOccluded1(rootScene, &shadowContext, &shadowRay);
#endif

                const float visibility = shadowRay.tfar < 0.0f ? 0.0f : 1.0f;
                pixelColor += albedo * visibility;
            }

            pixelColor = pixelColor * (1.0f / static_cast<float>(sampleCount));
            const size_t idx = (static_cast<size_t>(y) * static_cast<size_t>(width) +
                                static_cast<size_t>(x)) *
                               4u;
            image[idx + 0] = static_cast<uint8_t>(Clamp01(pixelColor.x) * 255.0f + 0.5f);
            image[idx + 1] = static_cast<uint8_t>(Clamp01(pixelColor.y) * 255.0f + 0.5f);
            image[idx + 2] = static_cast<uint8_t>(Clamp01(pixelColor.z) * 255.0f + 0.5f);
            image[idx + 3] = 255u;
        }
    }
    return image;
}
} // namespace

int main(int argc, char **argv)
{
#if !defined(WITH_EMBREE)
    (void)argc;
    (void)argv;
    fprintf(stderr, "Embree support is required for yeoubi_embree_harness.\n");
    return 1;
#else
    const CliOptions options = ParseCli(argc, argv);

    ScenePool scenePool;
    LoadUSDScene(&scenePool, options.inputPath);
    if (scenePool.scenes.empty() || scenePool.rootSceneIndex >= scenePool.scenes.size())
    {
        fprintf(stderr, "Failed to load USD scene: %s\n", options.inputPath.c_str());
        return 1;
    }

    ScenePool flattened;
    std::string flattenError;
    if (!FlattenScenePoolToRootChildren(&scenePool, &flattened, &flattenError))
    {
        fprintf(stderr, "Failed to flatten USD ScenePool: %s\n", flattenError.c_str());
        return 1;
    }
    if (flattened.scenes.empty() || flattened.rootSceneIndex >= flattened.scenes.size())
    {
        fprintf(stderr, "Flattened ScenePool invalid.\n");
        return 1;
    }

    std::vector<SceneMeshUploadRef> meshUploadRefs;
    CollectScenePoolMeshUploadRefs(&flattened, &meshUploadRefs);
    std::vector<float3> albedoByRefIndex(meshUploadRefs.size(), make_float3(0.7f));
    for (uint32_t i = 0; i < albedoByRefIndex.size(); i++)
    {
        albedoByRefIndex[i] = MakeFlatAlbedo(i);
    }

    CPUDevice cpuDevice;
    for (const std::unique_ptr<Scene> &scenePtr : flattened.scenes)
    {
        cpuDevice.BuildBVH(scenePtr.get());
    }

    Scene *rootScene = flattened.scenes[flattened.rootSceneIndex].get();
    RTCScene embreeScene = (RTCScene)rootScene->bvhHandle;
    if (!embreeScene)
    {
        fprintf(stderr, "Root Embree scene missing.\n");
        return 1;
    }

    const std::optional<RenderCamera> renderCamera = BuildUsdRenderCamera(flattened.camera);
    if (!renderCamera.has_value())
    {
        fprintf(stderr, "USD camera missing/invalid: %s\n", options.inputPath.c_str());
        return 1;
    }

    Bounds sceneBounds = {};
    ComputeWorldBounds(flattened, sceneBounds);
    const float diagonal = sceneBounds.valid ? length(sceneBounds.max - sceneBounds.min) : 1.0f;
    const float aoBias = std::max(1e-4f, diagonal * 0.002f);
    const float aoDistance =
        options.aoMaxDistance > 0.0f ? options.aoMaxDistance : std::max(0.1f, diagonal * 0.25f);

    printf("embree_harness: render start (%dx%d, spp=%d, integrator=%s)\n",
           renderCamera->width,
           renderCamera->height,
           options.spp,
           options.integrator == IntegratorType::AO ? "ao" : "primary");
    std::vector<uint8_t> image = RenderImage(embreeScene,
                                             *renderCamera,
                                             options.integrator,
                                             options.spp,
                                             aoBias,
                                             aoDistance,
                                             albedoByRefIndex);
    if (!SavePNG(options.outputPath.c_str(), image, renderCamera->width, renderCamera->height))
    {
        fprintf(stderr, "Failed to write PNG: %s\n", options.outputPath.c_str());
        return 1;
    }
    printf("embree_harness: wrote %s\n", options.outputPath.c_str());
    return 0;
#endif
}
