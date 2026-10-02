#pragma once

#include "core/assert.h"
#include "kernel/launch/shading_context.h"
#include "kernel/shading/shading.h"
#include "kernel/shading/shading_helpers.h"
#include "kernel/geometry/cluster.h"
#include "kernel/scene.h"

namespace ybi
{

YBI_GPU bool FindAttribute(const ShadingContext &context,
                           AttributeId id,
                           uint32_t offset,
                           uint32_t count,
                           int &attributeIndex)
{
    for (uint32_t index = offset; index < offset + count; index++)
    {
        const DeviceAttribute &deviceAttribute = context.scene.deviceAttributes[index];
        if (deviceAttribute.id == id)
        {
            attributeIndex = index;
            return true;
        }
    }
    return false;
}

YBI_GPU void ShaderAttributeValueVector2(const ShadingContext &context,
                                         const ShadingInfo &info,
                                         float *registers,
                                         uint32_t &nodeOffset)
{
    (void)registers;

    AttributeId id = static_cast<AttributeId>(context.scene.shaderNodes[nodeOffset++]);
    uint32_t outputSlot = context.scene.shaderNodes[nodeOffset++];

    const Object &object = context.scene.objects[info.objectId];

    int attributeIndex;
    bool found =
        FindAttribute(context, id, object.primvarOffset, object.primvarCount, attributeIndex);

    YBI_ASSERT(found && "Requested material attribute was not found");

    const DeviceAttribute &attribute = context.scene.deviceAttributes[attributeIndex];
    Vec2f values[3];

    if (IsClusteredTriangle(object.clusterOffset, info.clusterId))
    {
        const uint32_t clusterIndex = object.clusterOffset + info.clusterId;
        const DeviceCluster &cluster = context.scene.compressedClusters[clusterIndex];
        const uint8_t *clusterData = context.scene.compressedClusterBytes + cluster.byteOffset;

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
