#pragma once

#include "core/assert.h"
#include "core/math.h"
#include "core/vec3f.h"
#include "render/bsdf/materialx/conductor_bsdf.h"
#include "render/bsdf/materialx/dielectric_bsdf.h"
#include "render/bsdf/materialx/oren_nayar_bsdf.h"
#include "render/bsdf/materialx/sheen.h"
#include "render/material_opcodes.h"
#include "scene/device_scene.h"
#include <cstdint>

namespace ybi
{

struct MaterialNodeSummary
{
    Vec3f albedo;
    Vec3f transmittance;
};

// The descriptors must be child-before-parent, with local child indices. The
// caller owns both the resolved parameter records and the per-node summaries.
// The root summary is summaries[nodeCount - 1]. viewDirection is already
// normalized and expressed in the material's local shading frame.
YBI_CPU_GPU bool CalculateTransmittance(const MaterialDescriptor *descriptors,
                                        uint32_t nodeCount,
                                        const BsdfClosure &closure,
                                        const Vec3f &viewDirection,
                                        MaterialNodeSummary *summaries,
                                        uint32_t summaryCapacity)
{
    YBI_ASSERT(descriptors && summaries && nodeCount > 0 && summaryCapacity >= nodeCount);

    const float nDotV = viewDirection.z;

    for (uint32_t nodeIndex = 0; nodeIndex < nodeCount; ++nodeIndex)
    {
        const MaterialDescriptor &descriptor = descriptors[nodeIndex];
        MaterialNodeSummary &summary = summaries[nodeIndex];
        if (descriptor.flags != 0)
            return false;

        switch (descriptor.materialType)
        {
            case MaterialOpcode::layer_bsdf:
            {
                if (descriptor.child0 >= nodeIndex || descriptor.child1 >= nodeIndex)
                    return false;
                const MaterialNodeSummary &top = summaries[descriptor.child0];
                const MaterialNodeSummary &base = summaries[descriptor.child1];
                summary.albedo = top.albedo + top.transmittance * base.albedo;
                summary.transmittance = top.transmittance * base.transmittance;
                break;
            }
            case MaterialOpcode::mix_bsdf:
            {
                if (descriptor.child0 >= nodeIndex || descriptor.child1 >= nodeIndex)
                    return false;
                float weight;
                if (!closure.LoadFloat(descriptor.parameterOffset, weight))
                    return false;
                const MaterialNodeSummary &bg = summaries[descriptor.child0];
                const MaterialNodeSummary &fg = summaries[descriptor.child1];
                summary.albedo = (1.f - weight) * bg.albedo + weight * fg.albedo;
                summary.transmittance =
                    (1.f - weight) * bg.transmittance + weight * fg.transmittance;
                break;
            }
            case MaterialOpcode::sheen_bsdf:
            {
                SheenBSDF sheen;
                if (!GetSheenParameters(descriptor, closure, sheen))
                    return false;
                const float physicalAlbedo = sheen.weight * SheenAlbedo(sheen.roughness, nDotV);
                summary.albedo = physicalAlbedo * Max(sheen.color, 0.f);
                summary.transmittance = Vec3f(1.f - physicalAlbedo);
                break;
            }
            case MaterialOpcode::dielectric_bsdf:
            {
                DielectricBsdf dielectric;
                if (!GetDielectricParameters(descriptor, closure, dielectric))
                    return false;
                const Vec3f physicalAlbedo = CalculateDielectricBsdfAlbedo(dielectric, nDotV);
                summary.albedo = physicalAlbedo * Max(dielectric.tint, 0.f);
                summary.transmittance = Vec3f(1.f) - physicalAlbedo;
                break;
            }
            case MaterialOpcode::oren_nayar_diffuse_bsdf:
            {
                OrenNayarBsdf diffuse;
                if (!GetOrenNayarParameters(descriptor, closure, diffuse))
                    return false;
                summary.albedo = diffuse.weight * Max(diffuse.color, 0.f) *
                                 OrenNayarDirectionalAlbedo(nDotV, diffuse.roughness);
                summary.transmittance = Vec3f(1.f - diffuse.weight);
                break;
            }
            case MaterialOpcode::conductor_bsdf:
            {
                ConductorBsdf conductor;
                if (!GetConductorParameters(descriptor, closure, conductor))
                    return false;
                summary.albedo = ConductorDirectionalAlbedo(conductor, nDotV);
                summary.transmittance = Vec3f(1.f - conductor.weight);
                break;
            }
            default:
                return false;
        }
    }
    return true;
}

// The sampling/evaluation prototype is intentionally deferred. Only the
// transmittance and directional-albedo summary pass is implemented here.
#if 0
YBI_CPU_GPU Vec3f EvaluateBsdf() {}

YBI_GPU void SelectLobe(const KernelParams &params,
                        const BsdfClosure &closure,
                        const Vec3f &viewDirection,
                        float u)
{
    uint32_t offset = 0;

    // Max Lobes
    float registers[32];

    uint32_t nodeCount = 0; // TODO
    for (uint32_t nodeIndex = 0; nodeIndex < nodeCount; nodeIndex++)
    {
        const MaterialDescriptor &materialDescriptor = params.scene.materialDescriptors[nodeIndex];

        switch (materialDescriptor.materialType)
        {
#if 0
            case MaterialOpcode::mix_bsdf:
            {
                // TODO: get somehow
                float mixWeight = 0.f;
                float weight0 = 0.f;
                float weight1 = 0.f;

                float reservoirWeight0 = weight0 * (1.f - mixWeight);
                float reservoirWeight1 = weight1 * mixWeight;

                float totalWeight = reservoirWeight0 + reservoirWeight1;
                float prob = reservoirWeight0 / totalWeight;

                // TODO: have to examine graph
                float transmittance0;
                float transmittance1;

                if (u < prob)
                {
                    u = u / prob;
                }
                else
                {
                    u = (u - prob) / (1.f - prob);
                }

                break;
            }
#endif
            case MaterialOpcode::layer_bsdf:
            {
                float directionalAlbedoTop;
                float directionalAlbedoBottom;
                Vec3f transmittanceTop;

                float probability =
                    directionalAlbedoTop /
                    (directionalAlbedoTop + directionalAlbedoBottom * transmittanceTop);

                if (u < probability)
                {
                    // Traverse to top layer
                    u = u / probability;
                }
                else
                {
                    // Traverse to bottom layer
                    u = (u - probability) / (1.f - probability);
                }
            }
            case MaterialOpcode::sheen_bsdf:
            {
                SheenBSDF sheen;
                GetSheenParameters(materialDescriptor, closure, sheen);
                float probability = SheenProbability(sheen, viewDirection);

                registers[sheen.probabilityIndex];
                float something = registers[sheen.probabilityIndex];

                break;
            }
        }
    }
}
#endif

} // namespace ybi
