#pragma once

#include "core/assert.h"
#include "device/kernel_params.h"
#include "render/shaders/shading.h"
#include "render/shaders/shading_helpers.h"
#include "scene/device_cluster.h"
#include "scene/device_scene.h"

namespace ybi
{

YBI_GPU bool FindAttribute(const KernelParams &params,
                           AttributeId id,
                           uint32_t offset,
                           uint32_t count,
                           int &attributeIndex)
{
    for (uint32_t index = offset; index < offset + count; index++)
    {
        const DeviceAttribute &deviceAttribute = params.scene.deviceAttributes[index];
        if (deviceAttribute.id == id)
        {
            attributeIndex = index;
            return true;
        }
    }
    return false;
}

YBI_GPU void ShaderAttributeValueVector2(const KernelParams &params,
                                         const ShadingInfo &info,
                                         float *registers,
                                         uint32_t &nodeOffset)
{
    (void)registers;

    AttributeId id = static_cast<AttributeId>(params.scene.shaderNodes[nodeOffset++]);
    uint32_t outputSlot = params.scene.shaderNodes[nodeOffset++];

    const Object &object = params.scene.objects[info.objectId];

    int attributeIndex;
    bool found =
        FindAttribute(params, id, object.primvarOffset, object.primvarCount, attributeIndex);

    YBI_ASSERT(found && "Requested material attribute was not found");

    const DeviceAttribute &attribute = params.scene.deviceAttributes[attributeIndex];
    Vec2f values[3];

    if (IsClusteredTriangle(object.clusterOffset, info.clusterId))
    {
        const uint32_t clusterIndex = object.clusterOffset + info.clusterId;
        const DeviceCluster &cluster = params.scene.compressedClusters[clusterIndex];
        const uint8_t *clusterData = params.scene.compressedClusterBytes + cluster.byteOffset;

        switch (attribute.attributeType)
        {
            case AttributeType::Vertex:
            {
                // TODO: this is redundant computation
                const Vec3u vertexIndices = cluster.DecodeTriangle(clusterData, info.primIndex);
                for (uint32_t i = 0; i < 3; i++)
                {
                    values[i] = cluster.DecodeAttribute<Vec2f>(
                        clusterData, attributeIndex - object.primvarOffset, vertexIndices[i]);
                }
                break;
            }
            case AttributeType::FaceVarying:
            {
                for (uint32_t i = 0; i < 3; i++)
                {
                    const uint32_t cornerIndex = 3u * info.primIndex + i;
                    values[i] = cluster.DecodeAttribute<Vec2f>(
                        clusterData, attributeIndex - object.primvarOffset, cornerIndex);
                }
                break;
            }
            default:
                YBI_ASSERT(0);
        }
    }
    else
    {
        // TODO: not supported yet
        YBI_ASSERT(0);
    }

    // TODO IMPORTANT: this value isn't actually being used because i'm defaulting to
    // mesh uvs in ShadingInfo

    Vec2f value = InterpolateTriangleValues(values, info.bary);
    StoreFloat2Register(registers, value, outputSlot);
}

} // namespace ybi
