#include "core/affine_transform.h"
#include "core/assert.h"
#include "core/frame.h"
#include "core/math.h"
#include "core/matrix4x4.h"
#include "core/octahedral.h"
#include "core/vec3f.h"
#include "kernel/backend/dispatch.h"
#include "kernel/launch/kernel_params.h"
#if defined(YBI_MATERIAL_BACKEND_OPENPBR) == defined(YBI_MATERIAL_BACKEND_MATERIALX)
#error Exactly one material backend must be selected for the OptiX kernel.
#endif
#if defined(YBI_MATERIAL_BACKEND_OPENPBR)
#include "kernel/bsdf/openpbr/api.h"
#include "kernel/bsdf/openpbr/init.h"
#endif
#include "kernel/backend/optix/intersection_programs.h"
#include "kernel/backend/optix/visibility_mask.h"
#include "kernel/geometry/attributes.h"
#include "kernel/geometry/clustered_triangle.h"
#include "kernel/geometry/triangle.h"
#include "kernel/integrator/mis.h"
#include "kernel/integrator/shade_hit.h"
#include "kernel/light/sample.h"
#include "kernel/light/visibility.h"
#include "kernel/raytracing/ray_cone.h"
#include "kernel/raytracing/ray_offset.h"
#include "kernel/raytracing/trace.h"
#include "kernel/sampling/mapping.h"
#include "kernel/sampling/rng.h"
#include "kernel/shading/image.h"
#include "kernel/shading/shading.h"
#include "kernel/state/types.h"
#include "kernel/texture/types.h"
#include "kernel/texture/virtual_sampling.h"
#include <ctime>
#include <optix_device.h>
#include <optix_types.h>
#include <vector_types.h>

namespace ybi
{

// #define YBI_DEBUG_LOD_COLOR 1

struct Ray
{
    Vec3f pos;
    Vec3f dir;
};

extern "C"
{
    __constant__ KernelParams params;
}

YBI_GPU static Vec3f GetLodDebugColor(float lod, uint32_t maxMip)
{
    const float normalizedLod = maxMip ? Clamp(lod / float(maxMip), 0.f, 1.f) : 0.f;

    if (normalizedLod < 0.25f)
    {
        return Lerp(Vec3f(0.f, 0.f, 1.f), Vec3f(0.f, 1.f, 1.f), normalizedLod * 4.f);
    }
    else if (normalizedLod < 0.5f)
    {
        return Lerp(Vec3f(0.f, 1.f, 1.f), Vec3f(0.f, 1.f, 0.f), (normalizedLod - 0.25f) * 4.f);
    }
    else if (normalizedLod < 0.75f)
    {
        return Lerp(Vec3f(0.f, 1.f, 0.f), Vec3f(1.f, 1.f, 0.f), (normalizedLod - 0.5f) * 4.f);
    }

    return Lerp(Vec3f(1.f, 1.f, 0.f), Vec3f(1.f, 0.f, 0.f), (normalizedLod - 0.75f) * 4.f);
}

#if 0
extern "C" __global__ void __raygen__kernel()
{
    const uint3 launchIndex = optixGetLaunchIndex();
    const uint3 launchDimensions = optixGetLaunchDimensions();
    const uint32_t rasterX = launchIndex.x + params.pixelOffsetX;
    const uint32_t rasterY = launchIndex.y + params.pixelOffsetY;

    RNG rng = RNG::Init(Vec2u(rasterX, rasterY), 0u);
    Vec3f pixel(float(rasterX) + 0.5f, float(rasterY) + 0.5f, 0.f);

    float filterU = rng.Uniform();
    float filterV = rng.Uniform();

    pixel.x += filterU - 0.5f;
    pixel.y += filterV - 0.5f;

    const int maxDepth = params.maxDepth;

    Matrix4x4 cameraFromRaster = params.cameraFromRaster;
    AffineTransform renderFromCamera = params.renderFromCamera;

    // Generate ray
    Vec4f cameraPos = cameraFromRaster * Vec4f(pixel, 1.f);
    cameraPos.xyz /= cameraPos.w;

    Vec3f rayOrigin = renderFromCamera.TransformPoint(Vec3f(0.f));
    Vec3f rayDirection = Normalize(renderFromCamera.TransformVector(cameraPos.xyz));
    float tMin = 0.f;
    float tMax = 1000000.f;
    float time = 0.f;
    uint32_t mask = 0xff;

    RayCone rayCone;
    rayCone.width = 0.f;
    rayCone.spreadAngle = atanf(2.f * tanf(0.5f * params.fovRadians) / float(params.imageHeight));

    OptixTraversableHandle bvhHandle = static_cast<OptixTraversableHandle>(params.bvhHandle);

    IntersectionInfo info;
    uint64_t ptr = reinterpret_cast<uint64_t>(&info);
    uint32_t low = ptr & 0xffffffffu;
    uint32_t high = ptr >> 32u;

    optixTrace(bvhHandle,
               make_float3(rayOrigin.x, rayOrigin.y, rayOrigin.z),
               make_float3(rayDirection.x, rayDirection.y, rayDirection.z),
               tMin,
               tMax,
               time,
               mask,
               OPTIX_RAY_FLAG_NONE,
               0,
               0,
               0,
               low,
               high);

    uchar4 color = make_uchar4(0, 0, 0, 255);
    if (info.type == INTERSECTION_TYPE_TRIANGLE)
    {
#if 1
        const Object &object = params.scene.objects[info.id];

        ShadingInfo shadingInfo;
        GetShadingInfo(object, info, shadingInfo);
        Vec3f *positions = shadingInfo.positions;
        Vec3f *normals = shadingInfo.normals;

        if (!shadingInfo.hasNormals)
        {
            Vec3f n = Normalize(Cross(positions[1] - positions[0], positions[2] - positions[0]));

            normals[0] = n;
            normals[1] = n;
            normals[2] = n;
        }

        Vec2f *uvs = shadingInfo.uvs;
        Vec2f uv = InterpolateTriangleValues(uvs, info.bary);

        if (shadingInfo.hasMaterial)
        {
            const MaterialHandle &handle = params.scene.materialHandles[shadingInfo.materialIndex];

#if defined(YBI_MATERIAL_BACKEND_OPENPBR)
            if (handle.type == MaterialHandleType::OpenPbr)
            {
                const DeviceMaterial2 &material = params.scene.deviceMaterial2s[handle.index];
                shadingInfo.TransformToWorld(info.objectToWorld);

                // Calculate LOD

                float surfaceSpreadAngle =
                    CalculatePrimaryHitUnifiedSurfaceSpreadAngle(rayCone,
                                                                 rayDirection,
                                                                 shadingInfo.shadingNormal,
                                                                 positions[0],
                                                                 positions[1],
                                                                 positions[2],
                                                                 normals[0],
                                                                 normals[1],
                                                                 normals[2]);
                PropagateRayCone(rayCone, surfaceSpreadAngle, info.t);
                shadingInfo.rayDirection = rayDirection;
                shadingInfo.rayCone = rayCone;

                const Vec3f baseColor =
                    EvaluateShaderVec3(params, shadingInfo, material.base_color);

                color =
                    make_uchar4(static_cast<unsigned char>(Clamp(baseColor.x, 0.f, 1.f) * 255.f),
                                static_cast<unsigned char>(Clamp(baseColor.y, 0.f, 1.f) * 255.f),
                                static_cast<unsigned char>(Clamp(baseColor.z, 0.f, 1.f) * 255.f),
                                255);
            }
            else
#endif
                // TODO: get rid of this legacy path
                if (handle.type == MaterialHandleType::Preview)
                {
                    DeviceMaterial &material = params.scene.materials[handle.index];
                    if (material.valid && shadingInfo.hasUvs)
                    {
                        shadingInfo.TransformToWorld(info.objectToWorld);

                        float surfaceSpreadAngle =
                            CalculatePrimaryHitUnifiedSurfaceSpreadAngle(rayCone,
                                                                         rayDirection,
                                                                         shadingInfo.shadingNormal,
                                                                         positions[0],
                                                                         positions[1],
                                                                         positions[2],
                                                                         normals[0],
                                                                         normals[1],
                                                                         normals[2]);
                        PropagateRayCone(rayCone, surfaceSpreadAngle, info.t);

                        shadingInfo.rayDirection = rayDirection;
                        shadingInfo.rayCone = rayCone;

                        const MaterialTextureHandle textureHandle = {material.udimBase,
                                                                     material.udimTableOffset};
                        const Vec4f value =
                            VirtualTextureSampleImage(params, textureHandle, shadingInfo);
                        color = make_uchar4(
                            static_cast<unsigned char>(Clamp(value.x, 0.f, 1.f) * 255.f),
                            static_cast<unsigned char>(Clamp(value.y, 0.f, 1.f) * 255.f),
                            static_cast<unsigned char>(Clamp(value.z, 0.f, 1.f) * 255.f),
                            static_cast<unsigned char>(Clamp(value.w, 0.f, 1.f) * 255.f));
                    }
#else
        uint32_t hash = info.clusterId;
        hash ^= hash >> 16u;
        hash *= 0x7feb352du;
        hash ^= hash >> 15u;
        hash *= 0x846ca68bu;
        hash ^= hash >> 16u;

        color.x = static_cast<unsigned char>(64u + (hash & 0xbfu));
        color.y = static_cast<unsigned char>(64u + ((hash >> 8u) & 0xbfu));
        color.z = static_cast<unsigned char>(64u + ((hash >> 16u) & 0xbfu));
#endif
                }
        }
    }

    uchar4 *image = reinterpret_cast<uchar4 *>(params.imageHandle);
    const uint32_t pixelIndex = launchIndex.x + launchIndex.y * launchDimensions.x;
    image[pixelIndex] = color;
}
#endif

YBI_GPU void RenderUnidirectional()
{
    const uint3 launchIndex = optixGetLaunchIndex();
    const uint3 launchDimensions = optixGetLaunchDimensions();
    const uint32_t rasterX = launchIndex.x + params.pixelOffsetX;
    const uint32_t rasterY = launchIndex.y + params.pixelOffsetY;
    const int maxDepth = Max(0, params.maxDepth);
    const OptixTraversableHandle bvhHandle = static_cast<OptixTraversableHandle>(params.bvhHandle);
    Vec3f sum(0.f);

    RNG rng = RNG::Init(Vec2u(rasterX, rasterY), params.sampleOffset);
    const float pixelX = float(rasterX) + rng.Uniform();
    const float pixelY = float(rasterY) + rng.Uniform();
    Vec4f cameraPos = params.cameraFromRaster * Vec4f(Vec3f(pixelX, pixelY, 0.f), 1.f);
    cameraPos.xyz /= cameraPos.w;
    Vec3f rayOrigin = params.renderFromCamera.TransformPoint(Vec3f(0.f));
    Vec3f rayDirection = Normalize(params.renderFromCamera.TransformVector(cameraPos.xyz));
    RayCone rayCone{};
    rayCone.spreadAngle = atanf(2.f * tanf(0.5f * params.fovRadians) / float(params.imageHeight));
    Vec3f radiance(0.f);
    Vec3f throughput(1.f);
    Vec3f previousSurfacePosition(0.f);
    float previousBsdfPdf = 0.f;
    uint32_t previousBsdfFlags = 0u;
    float etaScale = 1.f;
    bool hasPreviousBsdfSample = false;

    for (int depth = 0; depth <= maxDepth;)
    {
        const IntersectionInfo info =
            TraceRay(bvhHandle, rayOrigin, rayDirection, 0.f, 1.0e30f);

        if (info.type == INTERSECTION_TYPE_NONE)
        {
            Vec3f environment(0.f);
            const bool useMIS = hasPreviousBsdfSample && !(previousBsdfFlags & BSDF_DELTA);
            for (uint32_t lightIndex = 0; lightIndex < params.scene.numEnvmapLights; lightIndex++)
            {
                const float misWeight =
                    useMIS ? PowerHeuristic(previousBsdfPdf,
                                            PDFEnvmapLight(params.scene, lightIndex, rayDirection))
                           : 1.f;
                environment +=
                    EvaluateEnvmapLight(params.scene.envmapLights[lightIndex], rayDirection) *
                    misWeight;
            }
            radiance += throughput * environment;
            break;
        }

        const Object &object = params.scene.objects[info.id];
        if (object.lightIndex >= 0)
        {
            const DeviceLight &light = params.scene.deviceLights[object.lightIndex];

            // Rect lights emit from local -Z. OptiX triangles are two-sided, so a
            // back-face hit passes through without contributing.
            const Vec3f localRayDirection = light.localFromWorld.TransformVector(rayDirection);
            if (localRayDirection.z > 0.f)
            {
                float misWeight = 1.f;
                if (hasPreviousBsdfSample && !(previousBsdfFlags & BSDF_DELTA))
                {
                    const Vec3f lightHitPosition = rayOrigin + info.t * rayDirection;
                    const float lightPdf = PDFLight(params.scene,
                                                    static_cast<uint32_t>(object.lightIndex),
                                                    previousSurfacePosition,
                                                    lightHitPosition);
                    misWeight = PowerHeuristic(previousBsdfPdf, lightPdf);
                }
                radiance += throughput * light.radiance * misWeight;
            }

            const Vec3f lightHitPosition = rayOrigin + info.t * rayDirection;
            const Vec3f lightTangentX = info.objectToWorld.TransformVector(Vec3f(1.f, 0.f, 0.f));
            const Vec3f lightTangentY = info.objectToWorld.TransformVector(Vec3f(0.f, 1.f, 0.f));
            Vec3f lightNormal = Normalize(Cross(lightTangentX, lightTangentY));
            if (Dot(lightNormal, rayDirection) < 0.f)
            {
                lightNormal = -lightNormal;
            }
            rayOrigin = OffsetRayOrigin(lightHitPosition, lightNormal);
            continue;
        }

        if (depth == maxDepth)
            break;

        ShadingInfo shadingInfo;
        GetShadingInfo(params.scene, object, info, shadingInfo);
        shadingInfo.TransformToWorld(info.objectToWorld);

        Vec3f geometricNormal = shadingInfo.geometricNormal;
        if (Dot(geometricNormal, rayDirection) > 0.f)
            geometricNormal = -geometricNormal;

        // if (Dot(geometricNormal, rayDirection) > 0.f)
        // {
        //     shadingInfo.shadingNormal = -shadingInfo.shadingNormal;
        // }

        UpdateRayCone(rayCone, rayDirection, info.t, shadingInfo, depth);
        shadingInfo.rayDirection = rayDirection;
        shadingInfo.rayCone = rayCone;

#if defined(YBI_MATERIAL_BACKEND_OPENPBR)
        OpenPbrBsdf bsdf;
#endif

#if defined(YBI_MATERIAL_BACKEND_OPENPBR)
        const OrthonormalBasis tangentFrame(shadingInfo.dpdu, shadingInfo.shadingNormal);
#endif
        const Vec3f worldViewDirection = -rayDirection;
#if defined(YBI_MATERIAL_BACKEND_MATERIALX)
        const DeviceMaterial3 *materialXMaterial = nullptr;
        float closureData[MAX_CLOSURE_SIZE];
        ClosureInitializationInfo closureInitializationInfo;
        closureInitializationInfo.geometricNormal = shadingInfo.geometricNormal;
        closureInitializationInfo.shadingNormal = shadingInfo.shadingNormal;
        closureInitializationInfo.dpdu = shadingInfo.dpdu;
        closureInitializationInfo.closureData = closureData;
        closureInitializationInfo.uv = shadingInfo.uv;
        closureInitializationInfo.rayDirection = rayDirection;
        closureInitializationInfo.rayConeWidth = rayCone.width;
        closureInitializationInfo.triangleArea =
            Length(Cross(shadingInfo.positions[1] - shadingInfo.positions[0],
                         shadingInfo.positions[2] - shadingInfo.positions[0]));
        closureInitializationInfo.uvArea = Abs((shadingInfo.uvs[1].x - shadingInfo.uvs[0].x) *
                                                   (shadingInfo.uvs[2].y - shadingInfo.uvs[0].y) -
                                               (shadingInfo.uvs[1].y - shadingInfo.uvs[0].y) *
                                                   (shadingInfo.uvs[2].x - shadingInfo.uvs[0].x));
        const ClosureEvaluationInfo &closureEvaluationInfo = closureInitializationInfo;
#endif

        bool passedAlphaTest = true;
        bool hasBsdf = false;

        if (shadingInfo.hasMaterial)
        {
            const MaterialHandle &handle = params.scene.materialHandles[shadingInfo.materialIndex];
#if defined(YBI_MATERIAL_BACKEND_OPENPBR)
            if (handle.type == MaterialHandleType::OpenPbr)
            {
                const DeviceMaterial2 &material = params.scene.deviceMaterial2s[handle.index];

                passedAlphaTest = InitializeOpenPbrBsdf(params,
                                                        material,
                                                        shadingInfo,
                                                        bsdf,
                                                        shadingInfo.shadingNormal,
                                                        worldViewDirection,
                                                        tangentFrame,
                                                        throughput,
                                                        rng.Uniform());
                hasBsdf = true;
            }
#else
            if (handle.type == MaterialHandleType::MaterialX)
            {
                materialXMaterial = &params.scene.deviceMaterial3s[handle.index];
                DeviceDispatch<void>(materialXMaterial->initializeCallableIndex,
                                     &params,
                                     &closureInitializationInfo);
                hasBsdf = true;
            }
#endif
        }

        if (!shadingInfo.hasMaterial || !hasBsdf || !passedAlphaTest)
        {
            rayOrigin = OffsetRayOrigin(shadingInfo.position, -geometricNormal);
            ++depth;
            continue;
        }

        // Sample Light
        if (GetLightCount(params.scene) > 0)
        {
            LightSample lightSample;
            if (SampleLight(params.scene,
                            shadingInfo.position,
                            Vec3f(rng.Uniform(), rng.Uniform2D()),
                            lightSample))
            {
                Vec3f sampledLightBsdfWeight;
                float sampledLightBsdfPdf;
#if defined(YBI_MATERIAL_BACKEND_MATERIALX)
                sampledLightBsdfWeight =
                    DeviceDispatch<Vec3f>(materialXMaterial->evaluateCallableIndex,
                                          &params,
                                          &closureEvaluationInfo,
                                          worldViewDirection,
                                          lightSample.direction);
                sampledLightBsdfPdf = DeviceDispatch<float>(materialXMaterial->pdfCallableIndex,
                                                            &params,
                                                            &closureEvaluationInfo,
                                                            worldViewDirection,
                                                            lightSample.direction);
#else
                sampledLightBsdfWeight =
                    EvaluateOpenPBR(bsdf, worldViewDirection, lightSample.direction);
                sampledLightBsdfPdf = PDFOpenPBR(bsdf, worldViewDirection, lightSample.direction);
#endif
                if (sampledLightBsdfWeight.MaxComponent() != 0.f)
                {
                    if (IsLightSampleVisible(
                            bvhHandle, lightSample, shadingInfo.position, geometricNormal))
                    {
                        const float misWeight = PowerHeuristic(lightSample.pdf, sampledLightBsdfPdf);
                        radiance += throughput * lightSample.radiance * sampledLightBsdfWeight *
                                    misWeight / lightSample.pdf;
                    }
                }
            }
        }

        Vec3f worldLightDirection, weight;
        bool success = false;
#if defined(YBI_MATERIAL_BACKEND_MATERIALX)
        BsdfSample bsdfSample;
        if (materialXMaterial)
        {
            success = DeviceDispatch<bool>(materialXMaterial->sampleCallableIndex,
                                           &params,
                                           &closureEvaluationInfo,
                                           worldViewDirection,
                                           &bsdfSample,
                                           Vec3f(rng.Uniform(), rng.Uniform2D()));

            if (success)
            {
                worldLightDirection = bsdfSample.direction;
                weight = bsdfSample.weight;
                etaScale *= Sqr(bsdfSample.eta);
                previousBsdfFlags = bsdfSample.flags;
            }
        }
#else
        success = SampleOpenPBR(
            bsdf, worldViewDirection, throughput, rng.Uniform2D(), worldLightDirection, weight);
#endif

        if (!success)
            break;

        throughput *= weight;

        // Russian roulette
        const float survivalProbability = (throughput * etaScale).MaxComponent();
        if (survivalProbability < 1.f && depth > 1)
        {
            if (rng.Uniform() >= survivalProbability)
            {
                break;
            }
            throughput /= survivalProbability;
        }

        const Vec3f nextDirection = worldLightDirection;

        previousSurfacePosition = shadingInfo.position;
#if defined(YBI_MATERIAL_BACKEND_MATERIALX)
        if (materialXMaterial)
        {
            previousBsdfPdf = bsdfSample.pdf;
        }
#else
        previousBsdfPdf = PDFOpenPBR(bsdf, worldViewDirection, worldLightDirection);
#endif
        hasPreviousBsdfSample = true;

        // Both vectors must be in world space when choosing the offset side.
        if (Dot(nextDirection, geometricNormal) < 0.f)
        {
            geometricNormal = -geometricNormal;
        }
        rayOrigin = OffsetRayOrigin(shadingInfo.position, geometricNormal);
        rayDirection = nextDirection;
        ++depth;
    }
    sum += radiance;

    Vec3f *image = params.kernelOutput;
    const uint32_t pixelIndex = launchIndex.x + launchIndex.y * launchDimensions.x;
    image[pixelIndex] = sum;
}

extern "C" __global__ void __raygen__unidirectional()
{
    RenderUnidirectional();
}

} // namespace ybi
