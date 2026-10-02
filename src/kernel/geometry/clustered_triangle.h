#pragma once

#include "core/assert.h"
#include "device/cuda_compat.h"
#include "kernel/shading/shading.h"
#include "kernel/state/types.h"
#include "kernel/scene.h"

namespace ybi
{

YBI_GPU bool GetClusteredTriangleNormals(const Object &object,
                                         const IntersectionInfo &info,
                                         const DeviceCluster &cluster,
                                         const uint8_t *clusterData,
                                         const Vec3u &vertexIndices,
                                         Vec3f normals[3])
{
    switch (object.normalType)
    {
        case AttributeType::FaceVarying:
        {
            for (int i = 0; i < 3; i++)
            {
                const uint32_t normalIndex = 3u * info.primIndex + static_cast<uint32_t>(i);
                normals[i] = cluster.DecodeNormal(clusterData, normalIndex);
            }
            return true;
        }
        case AttributeType::Vertex:
        {
            for (uint32_t i = 0; i < 3u; i++)
            {
                normals[i] = cluster.DecodeNormal(clusterData, vertexIndices[i]);
            }
            return true;
        }
        case AttributeType::None:
        {
            return false;
        }
        default:
            YBI_ASSERT(0);
            return false;
    }
}

YBI_GPU void GetClusteredTrianglePositions(const DeviceCluster &cluster,
                                           const uint8_t *clusterData,
                                           const Vec3u &vertexIndices,
                                           Vec3f positions[3])
{
    for (uint32_t i = 0; i < 3u; i++)
    {
        positions[i] = cluster.DecodePosition(clusterData, vertexIndices[i]);
    }
}

YBI_GPU bool GetClusteredTriangleUVs(const Object &object,
                                     const IntersectionInfo &info,
                                     const DeviceCluster &cluster,
                                     const uint8_t *clusterData,
                                     const Vec3u &vertexIndices,
                                     Vec2f uvs[3])
{
    switch (object.uvType)
    {
        case AttributeType::FaceVarying:
        {
            for (uint32_t i = 0; i < 3u; i++)
            {
                const uint32_t uvIndex = 3u * info.primIndex + i;
                uvs[i] = cluster.DecodeUV(clusterData, uvIndex);
            }
            break;
        }
        case AttributeType::Vertex:
        {
            for (uint32_t i = 0; i < 3u; i++)
            {
                uvs[i] = cluster.DecodeUV(clusterData, vertexIndices[i]);
            }
            break;
        }
        case AttributeType::None:
        {
            return false;
        }
        default:
            YBI_ASSERT(0);
    }

    return true;
}

YBI_GPU bool GetClusteredTriangleMaterial(const IntersectionInfo &info,
                                          const DeviceCluster &cluster,
                                          const uint8_t *clusterData,
                                          uint32_t &materialIndex)
{
    materialIndex = cluster.DecodeMaterialID(clusterData, info.primIndex);
    if (materialIndex == DEVICE_CLUSTER_INVALID_MATERIAL_ID)
    {
        return false;
    }

    return true;
}

YBI_GPU void GetClusteredShadingInfo(const Object &object,
                                     const IntersectionInfo &info,
                                     const DeviceCluster &cluster,
                                     const uint8_t *clusterData,
                                     const Vec3u &vertexIndices,
                                     ShadingInfo &shadingInfo)
{
    GetClusteredTrianglePositions(cluster, clusterData, vertexIndices, shadingInfo.positions);
    shadingInfo.hasNormals = GetClusteredTriangleNormals(
        object, info, cluster, clusterData, vertexIndices, shadingInfo.normals);
    shadingInfo.hasUvs = GetClusteredTriangleUVs(
        object, info, cluster, clusterData, vertexIndices, shadingInfo.uvs);

    if (object.materialIndex != -1)
    {
        shadingInfo.materialIndex = object.materialIndex;
        shadingInfo.hasMaterial = true;
    }
    else
    {
        shadingInfo.hasMaterial =
            GetClusteredTriangleMaterial(info, cluster, clusterData, shadingInfo.materialIndex);
    }
}

} // namespace ybi
