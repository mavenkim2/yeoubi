#pragma once

#include "render/color_transform.h"
#include "render/integrator_bsdf.h"
#include "render/integrator_light.h"
#include "render/integrator_texture.h"
#include "render/shading_core.h"

#include <stdio.h>

namespace ybi
{
namespace render
{
namespace integrator
{

YBI_DEVICE const LaunchParams::InstanceGeomRef *GetGeomRefs(const LaunchParams &params)
{
    if (params.instanceGeomRefs == 0ull || params.instanceGeomRefCount <= 0)
    {
        return nullptr;
    }
    return reinterpret_cast<const LaunchParams::InstanceGeomRef *>(params.instanceGeomRefs);
}

YBI_DEVICE const PackedMaterial *GetPackedMaterials(const LaunchParams &params)
{
    if (params.materials == 0ull || params.materialCount <= 0)
    {
        return nullptr;
    }
    return reinterpret_cast<const PackedMaterial *>(params.materials);
}

YBI_DEVICE EvaluatedMaterial DefaultMaterial()
{
    return {};
}

YBI_DEVICE bool DebugTraceEnabled(const LaunchParams &params, const UInt2 &launchIndex)
{
    return params.singlePixelEnabled != 0 &&
           launchIndex.x == static_cast<unsigned int>(params.singlePixelX) &&
           launchIndex.y == static_cast<unsigned int>(params.singlePixelY);
}

template <typename State>
YBI_DEVICE bool DebugTraceEnabled(const State &state)
{
    return DebugTraceEnabled(state.Params(), state.LaunchIndex());
}

YBI_DEVICE void DebugPrintLightContribution(const LaunchParams &params,
                                            const UInt2 &launchIndex,
                                            const char *label,
                                            int depth,
                                            int lightIndex,
                                            unsigned int lightType,
                                            float lightPdf,
                                            float bsdfPdf,
                                            float misWeight,
                                            const Vec3 &throughput,
                                            const Vec3 &lightRadiance,
                                            const Vec3 &bsdfValue,
                                            const Vec3 &shadowTr,
                                            const Vec3 &contribution)
{
    if (!DebugTraceEnabled(params, launchIndex))
    {
        return;
    }

    printf("debug %s pixel=(%u,%u) spp=%d depth=%d light=%d type=%s lightPdf=%.9g bsdfPdf=%.9g mis=%.9g\n",
           label,
           launchIndex.x,
           launchIndex.y,
           params.currentSpp,
           depth,
           lightIndex,
           LightTypeName(lightType),
           lightPdf,
           bsdfPdf,
           misWeight);
    printf("  throughput=(%.6f %.6f %.6f) Li=(%.6f %.6f %.6f) f=(%.6f %.6f %.6f) Tr=(%.6f %.6f %.6f) contrib=(%.6f %.6f %.6f)\n",
           throughput.x,
           throughput.y,
           throughput.z,
           lightRadiance.x,
           lightRadiance.y,
           lightRadiance.z,
           bsdfValue.x,
           bsdfValue.y,
           bsdfValue.z,
           shadowTr.x,
           shadowTr.y,
           shadowTr.z,
           contribution.x,
           contribution.y,
           contribution.z);
}

YBI_DEVICE void DebugPrintEmissionContribution(const LaunchParams &params,
                                               const UInt2 &launchIndex,
                                               int depth,
                                               const Vec3 &throughput,
                                               const Vec3 &emission,
                                               const Vec3 &contribution)
{
    if (!DebugTraceEnabled(params, launchIndex))
    {
        return;
    }

    printf("debug emissive-hit pixel=(%u,%u) spp=%d depth=%d\n",
           launchIndex.x,
           launchIndex.y,
           params.currentSpp,
           depth);
    printf("  throughput=(%.6f %.6f %.6f) emission=(%.6f %.6f %.6f) contrib=(%.6f %.6f %.6f)\n",
           throughput.x,
           throughput.y,
           throughput.z,
           emission.x,
           emission.y,
           emission.z,
           contribution.x,
           contribution.y,
           contribution.z);
}

YBI_DEVICE void DebugPrintBsdfSample(const LaunchParams &params,
                                     const UInt2 &launchIndex,
                                     int depth,
                                     const BsdfSample &bsdf,
                                     const Vec3 &throughput)
{
    if (!DebugTraceEnabled(params, launchIndex))
    {
        return;
    }

    printf("debug bsdf pixel=(%u,%u) spp=%d depth=%d pdf=%.9g skipNeeMis=%d transmission=%d\n",
           launchIndex.x,
           launchIndex.y,
           params.currentSpp,
           depth,
           bsdf.pdf,
           bsdf.skipNeeMis ? 1 : 0,
           bsdf.transmission ? 1 : 0);
    printf("  wi=(%.6f %.6f %.6f) weight=(%.6f %.6f %.6f) throughputNext=(%.6f %.6f %.6f)\n",
           bsdf.wi.x,
           bsdf.wi.y,
           bsdf.wi.z,
           bsdf.weight.x,
           bsdf.weight.y,
           bsdf.weight.z,
           throughput.x,
           throughput.y,
           throughput.z);
}

YBI_DEVICE EvaluatedMaterial LoadEvaluatedMaterial(const LaunchParams &params, int materialIndex)
{
    EvaluatedMaterial material = DefaultMaterial();
    const PackedMaterial *materials = GetPackedMaterials(params);
    if (!materials || materialIndex < 0 || materialIndex >= params.materialCount)
    {
        return material;
    }

    const PackedMaterial &src = materials[materialIndex];
    material.baseColor = Vec3(src.baseColor.x, src.baseColor.y, src.baseColor.z);
    material.emissiveColor = Vec3(src.emissiveColor.x, src.emissiveColor.y, src.emissiveColor.z);
    material.specularColor = Vec3(src.specularColor.x, src.specularColor.y, src.specularColor.z);
    material.roughness = src.roughness;
    material.metallic = src.metallic;
    material.ior = src.ior;
    material.opacity = src.opacity;
    material.clearcoat = src.clearcoat;
    material.clearcoatRoughness = src.clearcoatRoughness;
    material.opacityThreshold = src.opacityThreshold;
    material.flags = src.flags;
    material.useSpecularWorkflow = src.useSpecularWorkflow;
    material.hasAuthoredUseSpecularWorkflow = src.hasAuthoredUseSpecularWorkflow;
    return material;
}

template <typename State>
YBI_DEVICE EvaluatedMaterial EvaluateMaterial(State &state,
                                              const LaunchParams::InstanceGeomRef &geomRef,
                                              const HitInfo &hit)
{
    EvaluatedMaterial material = LoadEvaluatedMaterial(state.Params(), geomRef.materialIndex);

    Vec4 sample = {};
    if (TrySampleMaterialTextureSemantic(state, geomRef, hit, kSemanticDiffuse, sample))
    {
        material.baseColor = Vec3(sample.x, sample.y, sample.z);
    }
    if (TrySampleMaterialTextureSemantic(state, geomRef, hit, kSemanticRoughness, sample))
    {
        material.roughness = sample.x;
    }
    if (TrySampleMaterialTextureSemantic(state, geomRef, hit, kSemanticMetallic, sample))
    {
        material.metallic = sample.x;
    }
    if (TrySampleMaterialTextureSemantic(state, geomRef, hit, kSemanticEmissive, sample))
    {
        material.emissiveColor = Vec3(sample.x, sample.y, sample.z);
    }
    if (TrySampleMaterialTextureSemantic(state, geomRef, hit, kSemanticSpecularColor, sample))
    {
        material.specularColor = Vec3(sample.x, sample.y, sample.z);
    }
    if (TrySampleMaterialTextureSemantic(state, geomRef, hit, kSemanticIor, sample))
    {
        material.ior = sample.x;
    }
    if (TrySampleMaterialTextureSemantic(state, geomRef, hit, kSemanticOpacity, sample))
    {
        material.opacity = sample.x;
    }
    if (TrySampleMaterialTextureSemantic(state, geomRef, hit, kSemanticClearcoat, sample))
    {
        material.clearcoat = sample.x;
    }
    if (TrySampleMaterialTextureSemantic(state, geomRef, hit, kSemanticClearcoatRoughness, sample))
    {
        material.clearcoatRoughness = sample.x;
    }

    material.baseColor = ClampVec3(material.baseColor, 0.0f, 1.0f);
    material.emissiveColor = ClampVec3(material.emissiveColor, 0.0f, 65504.0f);
    material.specularColor = ClampVec3(material.specularColor, 0.0f, 1.0f);
    material.roughness = MaxFloat(MinFloat(material.roughness, 1.0f), 0.0f);
    material.metallic = Clamp(material.metallic, 0.0f, 1.0f);
    material.ior = MaxFloat(material.ior, 1.0f);
    material.opacity = Clamp(material.opacity, 0.0f, 1.0f);
    material.clearcoat = Clamp(material.clearcoat, 0.0f, 1.0f);
    material.clearcoatRoughness = Clamp(material.clearcoatRoughness, 0.0f, 1.0f);
    material.opacityThreshold = Clamp(material.opacityThreshold, 0.0f, 1.0f);
    material.flags = 0u;
    if (material.emissiveColor.x > 0.0f || material.emissiveColor.y > 0.0f ||
        material.emissiveColor.z > 0.0f)
    {
        material.flags |= MATERIAL_FLAG_HAS_EMISSION;
    }
    return material;
}

template <typename State>
YBI_DEVICE void DebugPrintEvaluatedMaterial(const State &state,
                                            const LaunchParams::InstanceGeomRef &geomRef,
                                            const HitInfo &hit,
                                            const EvaluatedMaterial &material,
                                            int depth)
{
    const LaunchParams &params = state.Params();
    const UInt2 launchIndex = state.LaunchIndex();
    if (!DebugTraceEnabled(params, launchIndex))
    {
        return;
    }

    printf("debug material pixel=(%u,%u) spp=%d depth=%d instance=%d prim=%d material=%d t=%.6f\n",
           launchIndex.x,
           launchIndex.y,
           params.currentSpp,
           depth,
           hit.instanceId,
           hit.primitiveIndex,
           geomRef.materialIndex,
           hit.t);
    printf("  baseColor=(%.6f %.6f %.6f) emissiveColor=(%.6f %.6f %.6f) specularColor=(%.6f %.6f %.6f)\n",
           material.baseColor.x,
           material.baseColor.y,
           material.baseColor.z,
           material.emissiveColor.x,
           material.emissiveColor.y,
           material.emissiveColor.z,
           material.specularColor.x,
           material.specularColor.y,
           material.specularColor.z);
    printf("  roughness=%.6f metallic=%.6f ior=%.6f opacity=%.6f opacityThreshold=%.6f\n",
           material.roughness,
           material.metallic,
           material.ior,
           material.opacity,
           material.opacityThreshold);
    printf("  clearcoat=%.6f clearcoatRoughness=%.6f flags=%u useSpecularWorkflow=%u hasAuthoredUseSpecularWorkflow=%u\n",
           material.clearcoat,
           material.clearcoatRoughness,
           material.flags,
           material.useSpecularWorkflow,
           material.hasAuthoredUseSpecularWorkflow);
}

YBI_DEVICE LaunchParams::InstanceGeomRef ResolveGeomRef(
    const LaunchParams &params, const LaunchParams::InstanceGeomRef *geomRefs, int instanceId)
{
    LaunchParams::InstanceGeomRef geomRef = {};
    geomRef.materialIndex = -1;
    if (geomRefs && instanceId >= 0 && instanceId < params.instanceGeomRefCount)
    {
        geomRef = geomRefs[instanceId];
    }
    return geomRef;
}

YBI_DEVICE bool ShouldAlphaCutout(const EvaluatedMaterial &material)
{
    return material.opacityThreshold > 0.0f && material.opacity < material.opacityThreshold;
}

YBI_DEVICE Vec3 ShadowTransmittance(const EvaluatedMaterial &material)
{
    const float transmission = Clamp(1.0f - material.opacity, 0.0f, 1.0f);
    return Vec3(transmission, transmission, transmission);
}

template <typename State>
YBI_DEVICE Vec3 TraceShadowTransmittance(State &state,
                                         const LaunchParams::InstanceGeomRef *geomRefs,
                                         const Vec3 &origin,
                                         const Vec3 &direction,
                                         float tMin,
                                         float tMax,
                                         int lightIndex)
{
    const LaunchParams &params = state.Params();
    Vec3 transmittance = Vec3(1.0f, 1.0f, 1.0f);
    Vec3 currentOrigin = origin;
    float currentMax = tMax;
    constexpr int kMaxShadowHits = 64;

    for (int shadowHit = 0; shadowHit < kMaxShadowHits; ++shadowHit)
    {
        HitInfo hit = {};
        if (!state.TraceClosest(currentOrigin, direction, tMin, currentMax, &hit))
        {
            return transmittance;
        }

        const float advance = hit.t;
        if (advance <= tMin)
        {
            return Vec3(0.0f, 0.0f, 0.0f);
        }

        const Vec3 hitPoint = hit.rayOrigin + hit.rayDir * hit.t;
        currentMax -= advance;
        if (currentMax <= tMin)
        {
            return transmittance;
        }

        if (lightIndex >= 0 && IsRefExcludedFromLightShadow(params, lightIndex, hit.instanceId))
        {
            currentOrigin = OffsetRayOrigin(hitPoint, hit.geomNormal, direction);
            continue;
        }

        const LaunchParams::InstanceGeomRef geomRef =
            ResolveGeomRef(params, geomRefs, hit.instanceId);
        const EvaluatedMaterial material = EvaluateMaterial(state, geomRef, hit);
        if (ShouldAlphaCutout(material))
        {
            currentOrigin = OffsetRayOrigin(hitPoint, hit.geomNormal, direction);
            continue;
        }

        const Vec3 shadowTr = ShadowTransmittance(material);
        if (MaxComponent(shadowTr) <= 0.0f)
        {
            return Vec3(0.0f, 0.0f, 0.0f);
        }

        transmittance = transmittance * shadowTr;
        if (MaxComponent(transmittance) <= 1.0e-3f)
        {
            return Vec3(0.0f, 0.0f, 0.0f);
        }

        currentOrigin = OffsetRayOrigin(hitPoint, hit.geomNormal, direction);
    }

    return Vec3(0.0f, 0.0f, 0.0f);
}

YBI_DEVICE Vec3 AccumulateAnalyticLightRadiance(const LaunchParams &params,
                                                const UInt2 &launchIndex,
                                                int depth,
                                                const Vec3 &throughput,
                                                const Vec3 &rayOrigin,
                                                const Vec3 &rayDir,
                                                float rayBias,
                                                float tMax,
                                                bool hasBsdfContext,
                                                bool skipNeeMis,
                                                float bsdfPdf)
{
    Vec3 lightRadiance = Vec3(0.0f, 0.0f, 0.0f);
    float tMin = rayBias;
    const PackedLight *lights = GetPackedLights(params);

    // Analytic lights are visible along the segment but do not block it.
    while (tMin < tMax)
    {
        LightRayHit lightHit = {};
        if (!TraceAnalyticLight(params, rayOrigin, rayDir, tMin, tMax, &lightHit))
        {
            break;
        }

        float misWeight = 1.0f;
        if (hasBsdfContext && !skipNeeMis && !lightHit.isDeltaLight)
        {
            misWeight = bsdfPdf / MaxFloat(bsdfPdf + lightHit.pdf, 1.0e-6f);
        }
        const Vec3 weightedRadiance = lightHit.radiance * misWeight;
        lightRadiance = lightRadiance + weightedRadiance;
        if (lights && lightHit.lightIndex >= 0 && lightHit.lightIndex < params.lightCount)
        {
            DebugPrintLightContribution(params,
                                        launchIndex,
                                        "light-hit",
                                        depth,
                                        lightHit.lightIndex,
                                        lights[lightHit.lightIndex].type,
                                        lightHit.pdf,
                                        bsdfPdf,
                                        misWeight,
                                        throughput,
                                        lightHit.radiance,
                                        Vec3(1.0f, 1.0f, 1.0f),
                                        Vec3(1.0f, 1.0f, 1.0f),
                                        throughput * weightedRadiance);
        }
        tMin = MaxFloat(lightHit.distance + rayBias, tMin + rayBias);
    }

    return lightRadiance;
}

template <typename State>
YBI_DEVICE Vec3 IntegratorPathTrace(State &state, const Vec3 &origin, const Vec3 &direction)
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

        if (depth > 0)
        {
            const float sceneDistance = hasSceneHit ? hit.t : 1.0e20f;
            const Vec3 lightRadiance = AccumulateAnalyticLightRadiance(params,
                                                                       launchIndex,
                                                                       depth,
                                                                       throughput,
                                                                       rayOrigin,
                                                                       rayDir,
                                                                       rayBias,
                                                                       sceneDistance,
                                                                       currentRayHasBsdfContext,
                                                                       currentRaySkipNeeMis,
                                                                       currentRayBsdfPdf);
            if (MaxComponent(lightRadiance) > 0.0f)
            {
                radiance = radiance + throughput * lightRadiance;
            }
        }

        if (!hasSceneHit)
        {
            const Vec3 missRadiance = EvaluateMissInfiniteLightRadiance(state,
                                                                        rayDir,
                                                                        currentRayHasBsdfContext,
                                                                        currentRaySkipNeeMis,
                                                                        currentRayBsdfPdf);
            const Vec3 contribution = throughput * missRadiance;
            radiance = radiance + contribution;
            const PackedLight *lights = GetPackedLights(params);
            const int domeLightIndex = FindDomeLightIndex(params);
            float domePdf = 0.0f;
            if (lights && domeLightIndex >= 0 && domeLightIndex < params.lightCount &&
                currentRayHasBsdfContext)
            {
                domePdf = DomeDirectionPdf(params, lights[domeLightIndex], rayDir);
            }
            DebugPrintLightContribution(params,
                                        launchIndex,
                                        "miss",
                                        depth,
                                        domeLightIndex,
                                        static_cast<unsigned int>(LightType::Dome),
                                        domePdf,
                                        currentRayBsdfPdf,
                                        1.0f,
                                        throughput,
                                        missRadiance,
                                        Vec3(1.0f, 1.0f, 1.0f),
                                        Vec3(1.0f, 1.0f, 1.0f),
                                        contribution);
            break;
        }

        const Vec3 hitPoint = hit.rayOrigin + hit.rayDir * hit.t;

        const LaunchParams::InstanceGeomRef geomRef =
            ResolveGeomRef(params, geomRefs, hit.instanceId);

        const EvaluatedMaterial material = EvaluateMaterial(state, geomRef, hit);
        const ShadingFrame shadingFrame = EvaluateShadingFrame(state, geomRef, hit);
        const Vec3 wo = -rayDir;
        const PreparedMaterial preparedMaterial =
            PrepareMaterialBsdf(material, shadingFrame, throughput, wo);

        if (ShouldAlphaCutout(preparedMaterial))
        {
            rayOrigin = OffsetRayOrigin(hitPoint, preparedMaterial.geomNormal, rayDir);
            depth--;
            continue;
        }

        DebugPrintEvaluatedMaterial(state, geomRef, hit, material, depth);

        if ((preparedMaterial.flags & MATERIAL_FLAG_HAS_EMISSION) != 0u)
        {
            const Vec3 contribution = throughput * preparedMaterial.emissiveColor;
            radiance = radiance + contribution;
            DebugPrintEmissionContribution(
                params, launchIndex, depth, throughput, preparedMaterial.emissiveColor, contribution);
        }

        DirectLightSample lightSample = {};
        if (SampleDirectLight(params, hitPoint, rngState, &lightSample))
        {
            const Vec3 shadowOrigin =
                OffsetRayOrigin(hitPoint, preparedMaterial.geomNormal, lightSample.wi);
            const float shadowMax = lightSample.distance >= 1.0e19f
                                        ? 1.0e20f
                                        : MaxFloat(lightSample.distance - rayBias, rayBias);
            const Vec3 shadowTr = TraceShadowTransmittance(state,
                                                           geomRefs,
                                                           shadowOrigin,
                                                           lightSample.wi,
                                                           rayBias,
                                                           shadowMax,
                                                           lightSample.lightIndex);
            if (MaxComponent(shadowTr) > 0.0f)
            {
                float bsdfPdf = 0.0f;
                const Vec3 f = EvaluateBsdf(preparedMaterial, lightSample.wi, &bsdfPdf);
                if (MaxComponent(f) > 0.0f)
                {
                    const float misWeight =
                        lightSample.isDeltaLight
                            ? 1.0f
                            : lightSample.pdf / MaxFloat(lightSample.pdf + bsdfPdf, 1.0e-6f);
                    const Vec3 contribution = throughput * f * shadowTr * lightSample.radiance *
                                              (misWeight / MaxFloat(lightSample.pdf, 1.0e-6f));
                    radiance = radiance + contribution;
                    const PackedLight *lights = GetPackedLights(params);
                    const unsigned int lightType =
                        (lights && lightSample.lightIndex >= 0 && lightSample.lightIndex < params.lightCount)
                            ? lights[lightSample.lightIndex].type
                            : 0u;
                    DebugPrintLightContribution(params,
                                                launchIndex,
                                                "direct",
                                                depth,
                                                lightSample.lightIndex,
                                                lightType,
                                                lightSample.pdf,
                                                bsdfPdf,
                                                misWeight,
                                                throughput,
                                                lightSample.radiance,
                                                f,
                                                shadowTr,
                                                contribution);
                }
            }
        }

        DirectLightSample domeSample = {};
        if (SampleDomeLight(state, rngState, &domeSample))
        {
            const Vec3 shadowOrigin =
                OffsetRayOrigin(hitPoint, preparedMaterial.geomNormal, domeSample.wi);
            const Vec3 shadowTr = TraceShadowTransmittance(state,
                                                           geomRefs,
                                                           shadowOrigin,
                                                           domeSample.wi,
                                                           rayBias,
                                                           1.0e20f,
                                                           domeSample.lightIndex);
            if (MaxComponent(shadowTr) > 0.0f)
            {
                float bsdfPdf = 0.0f;
                const Vec3 f = EvaluateBsdf(preparedMaterial, domeSample.wi, &bsdfPdf);
                if (MaxComponent(f) > 0.0f)
                {
                    const float misWeight =
                        domeSample.pdf / MaxFloat(domeSample.pdf + bsdfPdf, 1.0e-6f);
                    const Vec3 contribution = throughput * f * shadowTr * domeSample.radiance *
                                              (misWeight / MaxFloat(domeSample.pdf, 1.0e-6f));
                    radiance = radiance + contribution;
                    const PackedLight *lights = GetPackedLights(params);
                    const unsigned int lightType =
                        (lights && domeSample.lightIndex >= 0 && domeSample.lightIndex < params.lightCount)
                            ? lights[domeSample.lightIndex].type
                            : 0u;
                    DebugPrintLightContribution(params,
                                                launchIndex,
                                                "dome",
                                                depth,
                                                domeSample.lightIndex,
                                                lightType,
                                                domeSample.pdf,
                                                bsdfPdf,
                                                misWeight,
                                                throughput,
                                                domeSample.radiance,
                                                f,
                                                shadowTr,
                                                contribution);
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
        DebugPrintBsdfSample(params, launchIndex, depth, bsdf, throughput);
        if (MaxComponent(throughput) <= 0.0f)
        {
            break;
        }

        rayOrigin = OffsetRayOrigin(hitPoint, preparedMaterial.geomNormal, bsdf.wi);
        rayDir = bsdf.wi;
        currentRayHasBsdfContext = true;
        currentRaySkipNeeMis = bsdf.skipNeeMis;
        currentRayBsdfPdf = bsdf.pdf;
    }

    return radiance;
}

} // namespace integrator
} // namespace render
} // namespace ybi
