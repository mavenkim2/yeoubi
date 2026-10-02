#pragma once

#include "core/assert.h"
#include "kernel/scene.h"

namespace ybi
{

YBI_GPU Vec3u GetTriangleVertexIndices(const DeviceScene &scene,
                                       const Object &object,
                                       const IntersectionInfo &info)
{
    const uint32_t offset = object.triangleIndexOffset + 3u * info.primIndex;
    return Vec3u(static_cast<uint32_t>(scene.triangleIndices[offset]),
                 static_cast<uint32_t>(scene.triangleIndices[offset + 1u]),
                 static_cast<uint32_t>(scene.triangleIndices[offset + 2u]));
}

YBI_GPU void GetTrianglePositions(const DeviceScene &scene,
                                  const Object &object,
                                  const Vec3u &vertexIndices,
                                  Vec3f positions[3])
{
    for (uint32_t i = 0; i < 3u; i++)
    {
        positions[i] = scene.positions[object.positionOffset + vertexIndices[i]];
    }
}

YBI_GPU bool GetTriangleNormals(const DeviceScene &scene,
                                const Object &object,
                                const IntersectionInfo &info,
                                const Vec3u &vertexIndices,
                                Vec3f normals[3])
{
    switch (object.normalType)
    {
        case AttributeType::FaceVarying:
        {
            for (int i = 0; i < 3; i++)
            {
                Vec3f objectNormal =
                    DecodeOctahedral(scene.normals[object.normalOffset + 3 * info.primIndex + i]);

                normals[i] = objectNormal;
            }
            return true;
        }
        case AttributeType::Vertex:
        {
            for (uint32_t i = 0; i < 3u; i++)
            {
                Vec3f objectNormal =
                    DecodeOctahedral(scene.normals[object.normalOffset + vertexIndices[i]]);
                normals[i] = objectNormal;
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

YBI_GPU bool GetTriangleUVs(const DeviceScene &scene,
                            const Object &object,
                            const IntersectionInfo &info,
                            const Vec3u &vertexIndices,
                            Vec2f uvs[3])
{
    switch (object.uvType)
    {
        case AttributeType::FaceVarying:
        {
            for (int i = 0; i < 3; i++)
            {
                const uint32_t cornerIndex = 3 * info.primIndex + i;
                int uvIndex = cornerIndex;
                if (object.uvIndexOffset != ~0u)
                {
                    uvIndex = scene.uvIndices[object.uvIndexOffset + cornerIndex];
                }

                uvs[i] = scene.uvs[object.uvOffset + uvIndex];
            }
            break;
        }
        case AttributeType::Vertex:
        {
            for (uint32_t i = 0; i < 3u; i++)
            {
                int uvIndex = static_cast<int>(vertexIndices[i]);
                if (object.uvIndexOffset != ~0u)
                {
                    uvIndex = scene.uvIndices[object.uvIndexOffset + vertexIndices[i]];
                }

                uvs[i] = scene.uvs[object.uvOffset + uvIndex];
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

YBI_GPU bool GetTriangleMaterial(const DeviceScene &scene,
                                 const Object &object,
                                 uint32_t triangleIndex,
                                 uint32_t &materialIndex)
{
    if (object.triangleMaterialOffset != ~0u)
    {
        materialIndex =
            scene.triangleMaterialIndices[object.triangleMaterialOffset + triangleIndex];
        return true;
    }
    else if (object.materialIndex != -1)
    {
        materialIndex = object.materialIndex;
        return true;
    }
    return false;
}

YBI_GPU void GetTriangleShadingInfo(const DeviceScene &scene,
                                    const Object &object,
                                    const IntersectionInfo &info,
                                    ShadingInfo &shadingInfo)
{
    const Vec3u vertexIndices = GetTriangleVertexIndices(scene, object, info);
    GetTrianglePositions(scene, object, vertexIndices, shadingInfo.positions);
    shadingInfo.hasNormals =
        GetTriangleNormals(scene, object, info, vertexIndices, shadingInfo.normals);
    shadingInfo.hasUvs = GetTriangleUVs(scene, object, info, vertexIndices, shadingInfo.uvs);
    shadingInfo.hasMaterial =
        GetTriangleMaterial(scene, object, info.primIndex, shadingInfo.materialIndex);
}

} // namespace ybi
