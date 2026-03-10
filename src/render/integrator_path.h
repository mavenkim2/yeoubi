#pragma once

#include "render/integrator_bsdf.h"
#include "render/integrator_light.h"
#include "render/integrator_texture.h"
#include "render/shading_core.h"

namespace ybi
{
namespace render
{
namespace integrator
{

YBI_INTEGRATOR_HD const LaunchParams::InstanceGeomRef *GetGeomRefs(const LaunchParams &params)
{
    if (params.instanceGeomRefs == 0ull || params.instanceGeomRefCount <= 0)
    {
        return nullptr;
    }
    return reinterpret_cast<const LaunchParams::InstanceGeomRef *>(params.instanceGeomRefs);
}

YBI_INTEGRATOR_HD const PackedMaterial *GetPackedMaterials(const LaunchParams &params)
{
    if (params.materials == 0ull || params.materialCount <= 0)
    {
        return nullptr;
    }
    return reinterpret_cast<const PackedMaterial *>(params.materials);
}

YBI_INTEGRATOR_HD EvaluatedMaterial DefaultMaterial()
{
    return {};
}

YBI_INTEGRATOR_HD float LinearToSrgb(float value)
{
    if (value <= 0.0f)
    {
        return 0.0f;
    }
    if (value <= 0.0031308f)
    {
        return value * 12.92f;
    }
    return 1.055f * powf(value, 1.0f / 2.4f) - 0.055f;
}

YBI_INTEGRATOR_HD Vec3 DisplayMapPathRadiance(const Vec3 &radiance)
{
    constexpr float kDisplayExposure = 16.0f;
    const Vec3 exposed = radiance * kDisplayExposure;
    const Vec3 mapped =
        Vec3(1.0f - expf(-exposed.x), 1.0f - expf(-exposed.y), 1.0f - expf(-exposed.z));
    return Vec3(LinearToSrgb(mapped.x), LinearToSrgb(mapped.y), LinearToSrgb(mapped.z));
}

YBI_INTEGRATOR_HD EvaluatedMaterial LoadEvaluatedMaterial(const LaunchParams &params,
                                                          int materialIndex)
{
    EvaluatedMaterial material = DefaultMaterial();
    const PackedMaterial *materials = GetPackedMaterials(params);
    if (!materials || materialIndex < 0 || materialIndex >= params.materialCount)
    {
        return material;
    }

    const PackedMaterial &src = materials[materialIndex];
    material.baseColor = Vec3(src.baseColor.x, src.baseColor.y, src.baseColor.z);
    material.emissiveColor =
        Vec3(src.emissiveColor.x, src.emissiveColor.y, src.emissiveColor.z);
    material.specularColor =
        Vec3(src.specularColor.x, src.specularColor.y, src.specularColor.z);
    material.roughness = src.roughness;
    material.metallic = src.metallic;
    material.ior = src.ior;
    material.opacity = src.opacity;
    material.clearcoat = src.clearcoat;
    material.clearcoatRoughness = src.clearcoatRoughness;
    material.opacityThreshold = src.opacityThreshold;
    material.flags = src.flags;
    material.useSpecularWorkflow = src.useSpecularWorkflow;
    return material;
}

template <typename State>
YBI_INTEGRATOR_HD EvaluatedMaterial EvaluateMaterial(State &state,
                                                     const LaunchParams::InstanceGeomRef &geomRef,
                                                     const HitInfo &hit)
{
    EvaluatedMaterial material = LoadEvaluatedMaterial(state.Params(), geomRef.materialIndex);

    Vec4 sample = {};
    if (TrySampleMaterialTextureSemantic(state, geomRef, hit, kSemanticDiffuse, sample))
    {
        material.baseColor = material.baseColor * Vec3(sample.x, sample.y, sample.z);
    }
    if (TrySampleMaterialTextureSemantic(state, geomRef, hit, kSemanticRoughness, sample))
    {
        material.roughness *= sample.x;
    }
    if (TrySampleMaterialTextureSemantic(state, geomRef, hit, kSemanticMetallic, sample))
    {
        material.metallic *= sample.x;
    }
    if (TrySampleMaterialTextureSemantic(state, geomRef, hit, kSemanticEmissive, sample))
    {
        material.emissiveColor = material.emissiveColor * Vec3(sample.x, sample.y, sample.z);
    }
    if (TrySampleMaterialTextureSemantic(state, geomRef, hit, kSemanticSpecularColor, sample))
    {
        material.specularColor = material.specularColor * Vec3(sample.x, sample.y, sample.z);
    }
    if (TrySampleMaterialTextureSemantic(state, geomRef, hit, kSemanticIor, sample))
    {
        material.ior *= MaxFloat(sample.x, 0.0f);
    }
    if (TrySampleMaterialTextureSemantic(state, geomRef, hit, kSemanticOpacity, sample))
    {
        material.opacity *= sample.x;
    }
    if (TrySampleMaterialTextureSemantic(state, geomRef, hit, kSemanticClearcoat, sample))
    {
        material.clearcoat *= sample.x;
    }
    if (TrySampleMaterialTextureSemantic(state, geomRef, hit, kSemanticClearcoatRoughness, sample))
    {
        material.clearcoatRoughness *= sample.x;
    }

    material.baseColor = ClampVec3(material.baseColor, 0.0f, 1.0f);
    material.emissiveColor = ClampVec3(material.emissiveColor, 0.0f, 65504.0f);
    material.specularColor = ClampVec3(material.specularColor, 0.0f, 1.0f);
    material.roughness = MaxFloat(MinFloat(material.roughness, 1.0f), 0.0f);
    material.metallic = Clamp01(material.metallic);
    material.ior = MaxFloat(material.ior, 1.0f);
    material.opacity = Clamp01(material.opacity);
    material.clearcoat = Clamp01(material.clearcoat);
    material.clearcoatRoughness = Clamp01(material.clearcoatRoughness);
    material.opacityThreshold = Clamp01(material.opacityThreshold);
    material.flags = 0u;
    if (material.emissiveColor.x > 0.0f || material.emissiveColor.y > 0.0f ||
        material.emissiveColor.z > 0.0f)
    {
        material.flags |= MATERIAL_FLAG_HAS_EMISSION;
    }
    return material;
}

template <typename State>
YBI_INTEGRATOR_HD uint32_t IntegratorPathTrace(State &state,
                                               const Vec3 &origin,
                                               const Vec3 &direction)
{
    const LaunchParams &params = state.Params();
    const LaunchParams::InstanceGeomRef *geomRefs = GetGeomRefs(params);
    Vec3 radiance = Vec3(0.0f, 0.0f, 0.0f);
    Vec3 throughput = Vec3(1.0f, 1.0f, 1.0f);
    Vec3 rayOrigin = origin;
    Vec3 rayDir = direction;
    const float rayBias = params.aoBias > 0.0f ? params.aoBias : 0.001f;
    bool currentRayHasBsdfContext = false;
    bool currentRaySkipNeeMis = false;
    float currentRayBsdfPdf = 0.0f;

    const UInt2 launchIndex = state.LaunchIndex();
    unsigned int rngState =
        Hash32((launchIndex.x + 1u) * 73856093u ^ (launchIndex.y + 1u) * 19349663u ^
               (static_cast<unsigned int>(params.currentSpp) + 1u) * 83492791u);

    for (int depth = 0; depth <= MaxInt(params.maxDepth, 0); ++depth)
    {
        HitInfo hit = {};
        const bool hasSceneHit = state.TraceClosest(rayOrigin, rayDir, rayBias, 1.0e20f, &hit);
        const float sceneDistance = hasSceneHit ? hit.t : 1.0e20f;

        LightRayHit lightHit = {};
        if (TraceAnalyticLight(params, rayOrigin, rayDir, sceneDistance, &lightHit))
        {
            float misWeight = 1.0f;
            if (currentRayHasBsdfContext && !currentRaySkipNeeMis && !lightHit.isDeltaLight)
            {
                misWeight =
                    currentRayBsdfPdf / MaxFloat(currentRayBsdfPdf + lightHit.pdf, 1.0e-6f);
            }
            radiance = radiance + throughput * (lightHit.radiance * misWeight);
            break;
        }

        if (!hasSceneHit)
        {
            float misWeight = 1.0f;
            if (currentRayHasBsdfContext && !currentRaySkipNeeMis &&
                FindDomeLightIndex(params) >= 0)
            {
                const float domePdf = DomeDirectionPdf();
                misWeight = currentRayBsdfPdf / MaxFloat(currentRayBsdfPdf + domePdf, 1.0e-6f);
            }
            radiance = radiance + throughput * (EvaluateEnvironmentRadiance(state, rayDir) * misWeight);
            break;
        }

        const Vec3 hitPoint = hit.rayOrigin + hit.rayDir * hit.t;

        LaunchParams::InstanceGeomRef geomRef = {};
        geomRef.materialIndex = -1;
        if (geomRefs && hit.instanceId >= 0 && hit.instanceId < params.instanceGeomRefCount)
        {
            geomRef = geomRefs[hit.instanceId];
        }

        const EvaluatedMaterial material = EvaluateMaterial(state, geomRef, hit);
        const ShadingFrame shadingFrame = EvaluateShadingFrame(state, geomRef, hit);
        const Vec3 wo = -rayDir;
        const PreparedMaterial preparedMaterial =
            PrepareMaterialBsdf(material, shadingFrame, throughput, wo);

        if (ShouldAlphaCutout(preparedMaterial))
        {
            rayOrigin = hitPoint + rayDir * rayBias;
            depth--;
            continue;
        }

        if ((preparedMaterial.flags & MATERIAL_FLAG_HAS_EMISSION) != 0u)
        {
            radiance = radiance + throughput * preparedMaterial.emissiveColor;
        }

        DirectLightSample lightSample = {};
        if (SampleDirectLight(params, hitPoint, rngState, &lightSample))
        {
            const float shadowNormalSign =
                Dot(preparedMaterial.geomNormal, lightSample.wi) >= 0.0f ? 1.0f : -1.0f;
            const Vec3 shadowOrigin =
                hitPoint + preparedMaterial.geomNormal * (rayBias * shadowNormalSign);
            const float shadowMax =
                lightSample.distance >= 1.0e19f
                    ? 1.0e20f
                    : MaxFloat(lightSample.distance - rayBias, rayBias);
            if (!state.TraceLightOcclusion(
                    shadowOrigin, lightSample.wi, rayBias, shadowMax, lightSample.lightIndex))
            {
                float bsdfPdf = 0.0f;
                const Vec3 f = EvaluateBsdf(preparedMaterial, lightSample.wi, &bsdfPdf);
                if (MaxComponent(f) > 0.0f)
                {
                    const float misWeight =
                        lightSample.isDeltaLight
                            ? 1.0f
                            : lightSample.pdf / MaxFloat(lightSample.pdf + bsdfPdf, 1.0e-6f);
                    radiance = radiance +
                               throughput * f * lightSample.radiance *
                                   (misWeight / MaxFloat(lightSample.pdf, 1.0e-6f));
                }
            }
        }

        DirectLightSample domeSample = {};
        if (SampleDomeLight(state, rngState, &domeSample))
        {
            const float shadowNormalSign =
                Dot(preparedMaterial.geomNormal, domeSample.wi) >= 0.0f ? 1.0f : -1.0f;
            const Vec3 shadowOrigin =
                hitPoint + preparedMaterial.geomNormal * (rayBias * shadowNormalSign);
            if (!state.TraceLightOcclusion(
                    shadowOrigin, domeSample.wi, rayBias, 1.0e20f, domeSample.lightIndex))
            {
                float bsdfPdf = 0.0f;
                const Vec3 f = EvaluateBsdf(preparedMaterial, domeSample.wi, &bsdfPdf);
                if (MaxComponent(f) > 0.0f)
                {
                    const float misWeight =
                        domeSample.pdf / MaxFloat(domeSample.pdf + bsdfPdf, 1.0e-6f);
                    radiance = radiance +
                               throughput * f * domeSample.radiance *
                                   (misWeight / MaxFloat(domeSample.pdf, 1.0e-6f));
                }
            }
        }

        if (depth >= MaxInt(params.maxDepth, 0))
        {
            break;
        }

        BsdfSample bsdf = {};
        if (!SampleBsdf(preparedMaterial, rngState, &bsdf))
        {
            break;
        }

        throughput = throughput * bsdf.weight;
        if (MaxComponent(throughput) <= 0.0f)
        {
            break;
        }

        const float normalSign = Dot(preparedMaterial.geomNormal, bsdf.wi) >= 0.0f ? 1.0f : -1.0f;
        rayOrigin = hitPoint + preparedMaterial.geomNormal * (rayBias * normalSign);
        rayDir = bsdf.wi;
        currentRayHasBsdfContext = true;
        currentRaySkipNeeMis = bsdf.skipNeeMis;
        currentRayBsdfPdf = bsdf.pdf;
    }

    const Vec3 display = DisplayMapPathRadiance(radiance);
    return ybi::render::PackRGB8(display.x, display.y, display.z);
}

} // namespace integrator
} // namespace render
} // namespace ybi
