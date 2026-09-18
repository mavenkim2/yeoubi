#pragma once

#include "device/kernel_params.h"
#include "render/bsdf/materialx/dielectric_bsdf.h"
#include "render/bsdf/materialx/sheen.h"
#include "render/material_opcodes.h"

namespace ybi
{

YBI_CPU_GPU void CalculateTransmittance(const KernelParams &params,
                                        const BsdfClosure &closure,
                                        const Vec3f &viewDirection,
                                        float *transmittances)
{
    uint32_t nodeCount = 5;

    // Nodes supported in reverse topological sort order
    for (int nodeIndex = nodeCount - 1; nodeIndex >= 0; nodeIndex--)
    {
        const MaterialDescriptor &materialDescriptor = params.scene.materialDescriptors[nodeIndex];

        switch (materialDescriptor.materialType)
        {
            case MaterialOpcode::layer_bsdf:
            {
                float transmittanceTop;
                float transmittanceBottom;

                float transmittance = transmittanceTop * transmittanceBottom;

                // TODO: store
                break;
            }
            case MaterialOpcode::mix_bsdf:
            {
                float transmittanceBg;
                float transmittanceFg;
                float weight;

                float transmittance = (1.f - weight) * transmittanceBg + weight * transmittanceFg;

                // TODO: store
                break;
            }
            case MaterialOpcode::sheen_bsdf:
            {
                SheenBSDF sheen;
                GetSheenParameters(materialDescriptor, closure, sheen);

                float transmittance = 1.f - SheenAlbedo(sheen, viewDirection);
                // TODO: store

                break;
            }
            case MaterialOpcode::dielectric_bsdf:
            {
                // TODO: transmitting dielectric bsdf
                float transmittance =
                    1.f - CalculateDielectricBsdfAlbedo(const DielectricBsdf &bsdf,
                                                        const int &viewDirection) break;
                break;
            }
        }
    }
}

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

} // namespace ybi
